# Small function to write to message and to log file.

Small function to write to message and to log file.

## Usage

``` r
write_message(message_str, log_file = NULL)
```

## Arguments

- message_str:

  A string to write as a message.

- log_file:

  A log filename.

## Value

A message and writes the message to the specified log file.

## Examples

``` r
write_message(message_str = "Finished Step 1", log_file = "log.file.txt")
#> Finished Step 1
```
