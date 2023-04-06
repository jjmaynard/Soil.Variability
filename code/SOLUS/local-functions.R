

getData <- function(bb, .variable, .SOLUS_variable, .rescale_value) {
  # init spatVect
  e <- vect(bb, crs = "OGC:CRS84")
  
  # download google cloud asset to temp file
  x <- sprintf("projects/ncss-30m-covariates/assets/SOLUS100m/%s", .SOLUS_variable) |> 
    gd_image_from_id() |> 
    gd_download(region = gd_bbox(ext(e))) |> 
    rast()
  
  # get gNATSGO WCS raster for extent
  m <- mukey.wcs(e, db = 'gNATSGO', quiet = TRUE)
  
  # extract RAT
  rat <- cats(m)[[1]]
  
  # thematic soil data via SDA
  p <-  get_SDA_property(property = .variable,
                         method = "Weighted Average", 
                         mukeys = as.integer(rat[[1]]),
                         miscellaneous_areas = FALSE,
                         include_minors = TRUE,
                         top_depth = 14,
                         bottom_depth = 16)
  
  # join aggregate soil data + RAT
  rat <- merge(rat, p, by = 'mukey', sort = FALSE, all.x = TRUE)
  levels(m) <- rat
  
  # convert mukey + RAT -> numeric raster
  activeCat(m) <- .variable
  
  # must specify the variable index
  # faster than catalyze
  m2 <- as.numeric(m, index = match(.variable, names(rat)))
  
  # check correct variable was transferred
  stopifnot(names(m2) == .variable)
  
  # crop gNATSGO to SOLUS extent
  m2 <- crop(m2, x)
  
  # resample to SOLUS grid
  # m4 <- resample(m3, x)
  
  # resample SOLUS to gNATSGO grid
  # only target layer
  solus <- resample(x[[3]], m2)
  
  # rescale to range used by gNATSGO as-needed
  solus <- solus * .rescale_value
  
  return(c(gNATSGO = m2, SOLUS = solus))
}

