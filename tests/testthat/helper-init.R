test_genbank_file <- function(filename = "NC_000932.gb"){
  
  file <- system.file(
    "extdata",
    filename,
    package = "cpDNABarcoder"
  )
  
  if(!nzchar(file)){
    
    stop(
      "Test GenBank file not found in package extdata: ",
      filename
    )
    
  }
  
  file
  
}