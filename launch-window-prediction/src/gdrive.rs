//! Автозавантаження карт рельєфу LOLA та SPICE/ANISE ядер з Google Drive
//! згідно з кореневим `drive-manifest.json` та локальним кешем `cache/gdrive/`.

use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::time::Duration;

use anyhow::{Context, Result, anyhow};
use serde::Deserialize;

use crate::terrain::GridMeta;

pub const KERNELS_LOGICAL_PREFIX: &str =
    "trajectory/Lunar-Window-Navigator/launch-window-prediction/data/kernels/";
pub const TERRAIN_LOGICAL_PREFIX: &str =
    "trajectory/Lunar-Window-Navigator/launch-window-prediction/data/terrain/";

#[derive(Debug, Clone, Deserialize)]
pub struct DriveManifest {
    #[serde(default)]
    pub files: BTreeMap<String, DriveFileEntry>,
}

#[derive(Debug, Clone, Deserialize)]
pub struct DriveFileEntry {
    #[serde(rename = "fileId")]
    pub file_id: Option<String>,
    pub size: Option<u64>,
    pub sha256: Option<String>,
}

pub fn prefer_local() -> bool {
    std::env::var("LWP_PREFER_LOCAL").ok().as_deref() == Some("1")
        || std::env::var("GDRIVE_PREFER_LOCAL").ok().as_deref() == Some("1")
}

pub fn find_manifest_path() -> Option<PathBuf> {
    for var in ["LWP_DRIVE_MANIFEST", "GDRIVE_MANIFEST_PATH"] {
        if let Ok(val) = std::env::var(var) {
            let trimmed = val.trim();
            if !trimmed.is_empty() {
                let p = PathBuf::from(trimmed);
                if p.is_file() {
                    return Some(p);
                }
            }
        }
    }

    let mut roots = Vec::new();
    if let Ok(cwd) = std::env::current_dir() {
        roots.push(cwd);
    }
    if let Ok(exe) = std::env::current_exe() {
        if let Some(dir) = exe.parent() {
            roots.push(dir.to_path_buf());
        }
    }

    for start in roots {
        let mut cur: Option<&Path> = Some(start.as_path());
        for _ in 0..8 {
            let Some(dir) = cur else { break };
            let candidate = dir.join("drive-manifest.json");
            if candidate.is_file() {
                return Some(candidate);
            }
            cur = dir.parent();
        }
    }
    None
}

pub fn load_manifest() -> Option<(PathBuf, DriveManifest)> {
    let path = find_manifest_path()?;
    let text = fs::read_to_string(&path).ok()?;
    let parsed: DriveManifest = serde_json::from_str(&text).ok()?;
    Some((path, parsed))
}

pub fn cache_root(manifest_path: &Path) -> PathBuf {
    for var in ["LWP_GDRIVE_CACHE_DIR", "GDRIVE_CACHE_DIR"] {
        if let Ok(val) = std::env::var(var) {
            let trimmed = val.trim();
            if !trimmed.is_empty() {
                return PathBuf::from(trimmed);
            }
        }
    }
    manifest_path
        .parent()
        .unwrap_or_else(|| Path::new("."))
        .join("cache")
        .join("gdrive")
}

pub fn has_kernels() -> bool {
    let Some((_, manifest)) = load_manifest() else {
        return false;
    };
    manifest.files.iter().any(|(k, v)| {
        k.starts_with(KERNELS_LOGICAL_PREFIX)
            && v.file_id.as_deref().is_some_and(|id| !id.trim().is_empty())
    })
}

pub fn has_terrain() -> bool {
    let Some((_, manifest)) = load_manifest() else {
        return false;
    };
    manifest.files.iter().any(|(k, v)| {
        k.starts_with(TERRAIN_LOGICAL_PREFIX)
            && v.file_id.as_deref().is_some_and(|id| !id.trim().is_empty())
    })
}

pub fn dir_has_kernel_files(dir: &Path) -> bool {
    let Ok(entries) = fs::read_dir(dir) else {
        return false;
    };
    const EXTS: [&str; 4] = ["bsp", "bpc", "pca", "epa"];
    entries.flatten().any(|e| {
        let p = e.path();
        p.is_file()
            && p.extension()
                .and_then(|s| s.to_str())
                .is_some_and(|ext| EXTS.contains(&ext.to_ascii_lowercase().as_str()))
    })
}

pub fn dir_has_terrain_images(dir: &Path) -> bool {
    let Ok(entries) = fs::read_dir(dir) else {
        return false;
    };
    entries.flatten().any(|e| {
        let p = e.path();
        p.is_file()
            && p.extension()
                .and_then(|s| s.to_str())
                .is_some_and(|ext| ext.eq_ignore_ascii_case("img"))
    })
}

fn extract_attr<'a>(fragment: &'a str, attr: &str) -> Option<&'a str> {
    let needle = format!("{attr}=\"");
    let start = fragment.find(&needle)? + needle.len();
    let rest = &fragment[start..];
    let end = rest.find('"')?;
    Some(&rest[..end])
}

fn decode_basic_entities(s: &str) -> String {
    s.replace("&amp;", "&")
        .replace("&quot;", "\"")
        .replace("&#39;", "'")
        .replace("&lt;", "<")
        .replace("&gt;", ">")
}

fn extract_confirm_url(html: &str, file_id: &str) -> Option<String> {
    if let Some(form_start) = html.find("id=\"download-form\"") {
        let form_slice = &html[form_start..];
        let action_search_start = html[..form_start].rfind("<form").unwrap_or(0);
        let action = extract_attr(&html[action_search_start..], "action")
            .map(decode_basic_entities)
            .unwrap_or_else(|| "https://drive.usercontent.google.com/download".to_string());

        let form_end = form_slice.find("</form>").unwrap_or(form_slice.len());
        let body = &form_slice[..form_end];
        let mut params = Vec::new();
        let mut cursor = 0;
        while let Some(rel) = body[cursor..].find("<input") {
            let tag_start = cursor + rel;
            let Some(rel_end) = body[tag_start..].find('>') else {
                break;
            };
            let tag = &body[tag_start..tag_start + rel_end + 1];
            if tag.contains("type=\"hidden\"") {
                if let (Some(name), Some(val)) = (extract_attr(tag, "name"), extract_attr(tag, "value")) {
                    params.push(format!("{}={}", name, decode_basic_entities(val)));
                }
            }
            cursor = tag_start + rel_end + 1;
        }
        if !params.is_empty() {
            return Some(format!("{action}?{}", params.join("&")));
        }
    }

    let confirm_pos = html.find("confirm=")?;
    let confirm_rest = &html[confirm_pos + "confirm=".len()..];
    let confirm: String = confirm_rest
        .chars()
        .take_while(|c| c.is_ascii_alphanumeric() || *c == '_' || *c == '-')
        .collect();
    if confirm.is_empty() {
        return None;
    }

    let uuid_extra = if let Some(uuid_pos) = html.find("uuid=") {
        let uuid_rest = &html[uuid_pos + "uuid=".len()..];
        let uuid: String = uuid_rest
            .chars()
            .take_while(|c| c.is_ascii_alphanumeric() || *c == '_' || *c == '-')
            .collect();
        if uuid.is_empty() {
            String::new()
        } else {
            format!("&uuid={uuid}")
        }
    } else {
        String::new()
    };

    Some(format!(
        "https://drive.usercontent.google.com/download?id={file_id}&export=download&confirm={confirm}{uuid_extra}"
    ))
}

fn is_already_cached(dest: &Path, expected_size: Option<u64>) -> bool {
    let Ok(meta) = fs::metadata(dest) else {
        return false;
    };
    if !meta.is_file() {
        return false;
    }
    match expected_size {
        Some(sz) => meta.len() == sz,
        None => meta.len() > 0,
    }
}

pub fn download_file(file_id: &str, dest: &Path, expected_size: Option<u64>, label: &str) -> Result<()> {
    if is_already_cached(dest, expected_size) {
        return Ok(());
    }

    if let Some(parent) = dest.parent() {
        fs::create_dir_all(parent)
            .with_context(|| format!("create cache dir {}", parent.display()))?;
    }

    let agent = ureq::Agent::new_with_config(
        ureq::config::Config::builder()
            .timeout_global(Some(Duration::from_secs(900)))
            .build(),
    );

    let api_key = std::env::var("GOOGLE_DRIVE_API_KEY")
        .ok()
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty());

    let primary_url = match api_key {
        Some(key) => format!(
            "https://www.googleapis.com/drive/v3/files/{file_id}?alt=media&key={key}"
        ),
        None => format!(
            "https://drive.usercontent.google.com/download?id={file_id}&export=download&confirm=t"
        ),
    };

    tracing::info!(file = label, file_id = file_id, "downloading asset from Google Drive");

    let mut resp = agent
        .get(&primary_url)
        .header("User-Agent", "NASA-Lunar-Navigator/1.0")
        .call()
        .with_context(|| format!("GET {primary_url}"))?;

    let ctype = resp
        .headers()
        .get("content-type")
        .and_then(|v| v.to_str().ok())
        .unwrap_or("")
        .to_ascii_lowercase();

    if ctype.contains("text/html") {
        let html = resp
            .into_body()
            .read_to_string()
            .context("read Drive confirmation HTML")?;
        let confirm_url = extract_confirm_url(&html, file_id).ok_or_else(|| {
            anyhow!(
                "Google Drive returned an HTML confirmation page for {label} ({file_id}). Ensure sharing is 'Anyone with the link'."
            )
        })?;
        resp = agent
            .get(&confirm_url)
            .header("User-Agent", "NASA-Lunar-Navigator/1.0")
            .call()
            .with_context(|| format!("GET confirm URL for {label}"))?;
        let retry_ctype = resp
            .headers()
            .get("content-type")
            .and_then(|v| v.to_str().ok())
            .unwrap_or("")
            .to_ascii_lowercase();
        if retry_ctype.contains("text/html") {
            return Err(anyhow!(
                "Google Drive blocked automated download for {label} ({file_id})"
            ));
        }
    }

    let tmp_path = dest.with_extension(format!(
        "{}.{}.tmp",
        dest.extension().and_then(|s| s.to_str()).unwrap_or("bin"),
        std::process::id()
    ));

    let write_res: Result<u64> = (|| {
        let mut reader = resp.into_body().into_reader();
        let mut out = File::create(&tmp_path)
            .with_context(|| format!("create temp file {}", tmp_path.display()))?;
        let mut buf = vec![0u8; 256 * 1024];
        let mut written = 0u64;
        loop {
            let n = reader.read(&mut buf).context("read Drive stream")?;
            if n == 0 {
                break;
            }
            out.write_all(&buf[..n]).context("write cache temp file")?;
            written += n as u64;
        }
        out.flush().context("flush cache temp file")?;
        Ok(written)
    })();

    let written = match write_res {
        Ok(n) => n,
        Err(e) => {
            let _ = fs::remove_file(&tmp_path);
            return Err(e);
        }
    };

    if let Some(exp) = expected_size {
        if written != exp {
            let _ = fs::remove_file(&tmp_path);
            return Err(anyhow!(
                "size mismatch downloading {label}: expected {exp} bytes, got {written}"
            ));
        }
    }

    if let Err(e) = fs::rename(&tmp_path, dest) {
        if is_already_cached(dest, expected_size) {
            let _ = fs::remove_file(&tmp_path);
        } else {
            let _ = fs::remove_file(&tmp_path);
            return Err(e).with_context(|| format!("rename {} -> {}", tmp_path.display(), dest.display()));
        }
    }

    Ok(())
}

pub fn ensure_kernels_dir() -> Result<PathBuf> {
    let (manifest_path, manifest) =
        load_manifest().ok_or_else(|| anyhow!("drive-manifest.json not found"))?;
    let target_dir = cache_root(&manifest_path)
        .join("trajectory")
        .join("Lunar-Window-Navigator")
        .join("launch-window-prediction")
        .join("data")
        .join("kernels");
    fs::create_dir_all(&target_dir)
        .with_context(|| format!("create {}", target_dir.display()))?;

    let mut count = 0usize;
    for (logical_path, entry) in &manifest.files {
        if !logical_path.starts_with(KERNELS_LOGICAL_PREFIX) {
            continue;
        }
        let Some(file_id) = entry.file_id.as_deref().map(str::trim).filter(|s| !s.is_empty()) else {
            continue;
        };
        let file_name = &logical_path[KERNELS_LOGICAL_PREFIX.len()..];
        let dest = target_dir.join(file_name);
        download_file(file_id, &dest, entry.size, logical_path)?;
        count += 1;
    }

    if count == 0 {
        return Err(anyhow!(
            "no kernel files configured in {}",
            manifest_path.display()
        ));
    }

    Ok(target_dir)
}

fn find_matching_label_text(lbl_dir: &Path, img_file_name: &str) -> Option<(PathBuf, String)> {
    let stem = Path::new(img_file_name).file_stem()?.to_str()?;
    let entries = fs::read_dir(lbl_dir).ok()?;
    for entry in entries.flatten() {
        let p = entry.path();
        if !p.is_file() {
            continue;
        }
        let is_lbl = p
            .extension()
            .and_then(|s| s.to_str())
            .is_some_and(|e| e.eq_ignore_ascii_case("lbl"));
        let same_stem = p
            .file_stem()
            .and_then(|s| s.to_str())
            .is_some_and(|s| s.eq_ignore_ascii_case(stem));
        if is_lbl && same_stem {
            if let Ok(text) = fs::read_to_string(&p) {
                return Some((p, text));
            }
        }
    }
    None
}

fn label_covers_site(label_text: &str, lat_deg: f64, lon_deg: f64) -> bool {
    let Ok(meta) = GridMeta::parse_label(label_text) else {
        return true;
    };
    if let crate::terrain::Projection::PolarStereographic { south, .. } = meta.projection {
        if south && lat_deg > 0.0 {
            return false;
        }
        if !south && lat_deg < 0.0 {
            return false;
        }
    }
    let Some((sample, line)) = meta.lat_lon_to_pixel(lat_deg, lon_deg) else {
        return false;
    };
    sample >= 0.0
        && sample < meta.samples as f64
        && line >= 0.0
        && line < meta.lines as f64
}

pub fn ensure_terrain_dir(site: Option<(f64, f64)>, default_lbl_dir: &Path) -> Result<PathBuf> {
    let (manifest_path, manifest) =
        load_manifest().ok_or_else(|| anyhow!("drive-manifest.json not found"))?;
    let target_dir = cache_root(&manifest_path)
        .join("trajectory")
        .join("Lunar-Window-Navigator")
        .join("launch-window-prediction")
        .join("data")
        .join("terrain");
    fs::create_dir_all(&target_dir)
        .with_context(|| format!("create {}", target_dir.display()))?;

    // Копіюємо .LBL файли метаданих із локальної папки проєкту в кеш
    if default_lbl_dir.is_dir() {
        for entry in fs::read_dir(default_lbl_dir)?.flatten() {
            let p = entry.path();
            if p.is_file()
                && p.extension()
                    .and_then(|s| s.to_str())
                    .is_some_and(|e| e.eq_ignore_ascii_case("lbl"))
            {
                if let Some(name) = p.file_name() {
                    let dst = target_dir.join(name);
                    let src_len = fs::metadata(&p).map(|m| m.len()).unwrap_or(0);
                    let dst_len = fs::metadata(&dst).map(|m| m.len()).unwrap_or(u64::MAX);
                    if !dst.exists() || src_len != dst_len {
                        fs::copy(&p, &dst).with_context(|| {
                            format!("copy label {} -> {}", p.display(), dst.display())
                        })?;
                    }
                }
            }
        }
    }

    let mut selected = 0usize;
    for (logical_path, entry) in &manifest.files {
        if !logical_path.starts_with(TERRAIN_LOGICAL_PREFIX) {
            continue;
        }
        let Some(file_id) = entry.file_id.as_deref().map(str::trim).filter(|s| !s.is_empty()) else {
            continue;
        };
        let file_name = &logical_path[TERRAIN_LOGICAL_PREFIX.len()..];
        let is_img = Path::new(file_name)
            .extension()
            .and_then(|s| s.to_str())
            .is_some_and(|e| e.eq_ignore_ascii_case("img"));

        if is_img {
            if let Some((lat, lon)) = site {
                if let Some((_, lbl_text)) = find_matching_label_text(&target_dir, file_name) {
                    if !label_covers_site(&lbl_text, lat, lon) {
                        tracing::info!(
                            map = file_name,
                            lat_deg = lat,
                            lon_deg = lon,
                            "skipping Google Drive terrain download (site outside map bounds)"
                        );
                        continue;
                    }
                }
            }
        }

        let dest = target_dir.join(file_name);
        download_file(file_id, &dest, entry.size, logical_path)?;
        selected += 1;
    }

    if selected == 0 {
        return Err(anyhow!(
            "no terrain maps downloaded from {}",
            manifest_path.display()
        ));
    }

    Ok(target_dir)
}
