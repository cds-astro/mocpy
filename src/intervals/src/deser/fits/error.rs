//! Module defining the MOC FITS read/write/parsing errors

use std::io::{self, BufRead};
use std::num::ParseIntError;

use quick_error::{quick_error, ResultExt};

// https://docs.rs/quick-error/2.0.1/quick_error/
quick_error! {
  #[derive(Debug)]
  pub enum FitsError {
    /// IO error
    ReadIo(context: String, err: io::Error) {
      context(context: &'a str, err: io::Error) -> (String::from(context), err)
      display("I/O error reading '{}': {}", context, err)
    }
    UnexpectedKeyword(context: String, expected: String, actual: String) {
      display("Wrong keyword in '{}'. Expected: '{}'. Actual: '{}'.", context, expected, actual)
    }
    ValueIndicatorNotFound(context: String, keyword_record: String) {
      display("Value indicator expected but not found in '{}' for keyword record '{}'.", context, keyword_record)
    }
    UnexpectedValue(context: String, keyword: String, expected: String, actual: String) {
      display("Wrong value in '{}' for keyword '{}'. Expected: '{}'. Actual: '{}'.", context, keyword, expected, actual)
    }
    UintValueNotFound(context: String, keyword_record: String) {
      display("Unsigned int value not found in '{}', keyword record '{}'.", context, keyword_record)
    }
    StringValueNotFound(context: String, keyword_record: String) {
      display("String value no found in '{}', keyword record '{}'", context, keyword_record)
    }
    WrongUintValue(context: String, err: ParseIntError) {
      context(context: String, err: ParseIntError) -> (context, err)
      display("Parse error in {}: {}", context, err)
    }
    MultipleKeyword(context: String, keyword: String) {
      display("FITS {} not valid. Keyword '{}' found more than once.", context, keyword)
    }
    MissingKeyword(context: String, keyword: String) {
      display("Missing keyword in {}. Not found: '{}'", context, keyword)
    }
    UncompatibleKeywordContent(context: String, keyword1: String, keyword2: String) {
      display("In {}, incompatible keyword values for {} and {}", context, keyword1, keyword2)
    }
  }
}
