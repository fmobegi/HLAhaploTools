use std::fmt;

/// Errors that stop a whole file from being read.
///
/// Problems confined to one sample (a bad typing, an unknown MAC code) are not
/// errors: the sample is skipped and recorded as an [`crate::Issue`] instead,
/// the same way Hapl-o-Mat discards a donor and carries on.
#[derive(Debug)]
pub enum Error {
    Io(std::io::Error),
    Format(String),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Io(e) => write!(f, "I/O error: {e}"),
            Error::Format(msg) => write!(f, "format error: {msg}"),
        }
    }
}

impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Error::Io(e) => Some(e),
            Error::Format(_) => None,
        }
    }
}

impl From<std::io::Error> for Error {
    fn from(e: std::io::Error) -> Self {
        Error::Io(e)
    }
}
