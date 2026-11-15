use std::path::PathBuf;

/// Единый тип ошибок библиотеки.
#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("I/O error at {}: {source}", path.display())]
    Io {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },

    #[error("TOML parse error in {}: {source}", path.display())]
    Toml {
        path: PathBuf,
        #[source]
        source: toml::de::Error,
    },

    #[error("invalid config: {0}")]
    Config(String),

    #[error("almanac error: {0}")]
    Almanac(String),
}

pub type Result<T> = std::result::Result<T, Error>;
