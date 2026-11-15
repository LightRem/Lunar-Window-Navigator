//! Загрузка эфемерид и ориентаций в `anise::Almanac`.

use std::path::{Path, PathBuf};

use anise::prelude::{Almanac, MetaAlmanac};

use crate::error::{Error, Result};

/// Откуда брать SPICE/ANISE-ядра.
#[derive(Debug, Clone)]
pub enum AlmanacSource {
    /// `MetaAlmanac::latest()`: скачивает и кэширует de440s, pck11.pca, moon_fk_de440.epa,
    /// moon_pa_de440_200625.bpc и earth_latest_high_prec.bpc (нужен доступ в интернет).
    Latest,
    /// Локальная папка с файлами `.bsp`, `.bpc`, `.pca`, `.epa` (например, `data/kernels`).
    LocalDir(PathBuf),
}

/// Загружает `Almanac`. Сам разбор файлов синхронный и блокирующий,
/// поэтому он вынесен из async-рантайма в `spawn_blocking`.
pub async fn load_almanac(source: AlmanacSource) -> Result<Almanac> {
    tokio::task::spawn_blocking(move || load_almanac_blocking(source))
        .await
        .map_err(|e| Error::Almanac(format!("loader task failed: {e}")))?
}

/// Синхронная версия [`load_almanac`] (для тестов и кода, который уже работает вне async).
pub fn load_almanac_blocking(source: AlmanacSource) -> Result<Almanac> {
    match source {
        AlmanacSource::Latest => MetaAlmanac::latest().map_err(|e| Error::Almanac(e.to_string())),
        AlmanacSource::LocalDir(dir) => load_local_dir(&dir),
    }
}

fn load_local_dir(dir: &Path) -> Result<Almanac> {
    const EXTENSIONS: [&str; 4] = ["bsp", "bpc", "pca", "epa"];

    let entries = std::fs::read_dir(dir)
        .map_err(|source| Error::Io { path: dir.to_path_buf(), source })?;

    let mut files: Vec<PathBuf> = entries
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| {
            p.extension()
                .and_then(|e| e.to_str())
                .is_some_and(|e| EXTENSIONS.contains(&e.to_ascii_lowercase().as_str()))
        })
        .collect();
    files.sort();

    if files.is_empty() {
        return Err(Error::Config(format!(
            "no kernel files ({}) found in {}",
            EXTENSIONS.join(", "),
            dir.display()
        )));
    }

    let mut almanac = Almanac::default();
    for file in files {
        let path = file
            .to_str()
            .ok_or_else(|| Error::Config(format!("non-UTF-8 path: {}", file.display())))?;
        almanac = almanac
            .load(path)
            .map_err(|e| Error::Almanac(format!("{}: {e}", file.display())))?;
    }
    Ok(almanac)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn empty_dir_is_an_error() {
        let dir = std::env::temp_dir().join("lwp-empty-kernels-test");
        std::fs::create_dir_all(&dir).unwrap();
        let err = load_almanac(AlmanacSource::LocalDir(dir)).await.err().unwrap();
        assert!(matches!(err, Error::Config(_)));
    }

    #[tokio::test]
    async fn missing_dir_is_an_io_error() {
        let dir = std::env::temp_dir().join("lwp-no-such-dir-xyz");
        let _ = std::fs::remove_dir_all(&dir);
        let err = load_almanac(AlmanacSource::LocalDir(dir)).await.err().unwrap();
        assert!(matches!(err, Error::Io { .. }));
    }
}
