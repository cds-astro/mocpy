//! Module defining the MOC FITS read/write/parsing errors

use std::io;
use std::num::ParseIntError;

use quick_error::quick_error;

// https://docs.rs/quick-error/2.0.1/quick_error/
quick_error! {
  #[derive(Debug)]
  pub enum FitsError {
    /// IO error
    Io(err: io::Error) {
      from()
      display("I/O error: {}", err)
    }
    UnexpectedKeyword(expected: String, actual: String) {
      display("Wrong keyword. Expected: '{}'. Actual: '{}'.", expected, actual)
    }
    ValueIndicatorNotFound(keyword_record: String) {
      display("Value indicator not found in keyword record '{}'.", keyword_record)
    }
    UnexpectedValue(keyword: String, expected: String, actual: String) {
      display("Wrong value for keyword '{}'. Expected: '{}'. Actual: '{}'.", keyword, expected, actual)
    }
    UintValueNotFound(keyword_record: String) {
      display("Unsigned int value not found in keyword record '{}'.", keyword_record)
    }
    StringValueNotFound(keyword_record: String) {
      display("String value no found in keyword record '{}'", keyword_record)
    }
    WrongUintValue(context: String, err: ParseIntError) {
      context(context: String, err: ParseIntError) -> (context, err)
      display("Parse {} error: {}", context, err)
    }
    MultipleKeyword(keyword: String) {
      display("FITS not valid. Multiple Keyword '{}'", keyword)
    }
    MissingKeyword(keyword: String) {
      display("Missing keyword '{}'", keyword)
    }
    UncompatibleKeywordContent(keyword1: String, keyword2: String) {
      display("Incompatible keyword values for {} and {}", keyword1, keyword2)
    }
    RemainingData {
      display("More data than the expected!")
    }
    PrematureEndOfData {
      display("Less data than expected!")
    }
    UnexpectedWrittenSize {
      display("Unexpected number of data written!")
    }
    UnexpectedDepth(depth: u8, depth_max: u8) {
      display("unexpected depth. Max expected: {}. Actual: {}", depth_max, depth)
    }
  }
}
