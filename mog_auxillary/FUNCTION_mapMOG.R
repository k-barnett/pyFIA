#this script provides a function to assess if forest is classified as mature/old growth per the definitions set forth by the US Forest Service using the US Forest Service's Forest Inventory and Analysis data
#these definitions are loosely defined in Pelz, et al. (2023) and Woodall, et al. (2023)
#these definitions are explicitly defined in "Mature and Old-Growth Forests: Definition, Identification, and Initial Inventory on Lands Managed by the Forest Service and Bureau of Land Management. Fulfillment of Executive Order 14072, Section 2(b)". Throughout this script, I refer to this document as "the MOG document".
mapMOG <- function(locale, source.path, interpolate = TRUE, recent = TRUE, API = "rFIA", resolution = 1590, power = 2){ #begin function
  
  #warnings and errors
  if(!(class(locale)[1] %in% "sf")){stop("locale must be an 'sf' object.")}
  if(is.null(source.path)){stop("source.path argument must be provided. Please provide the file path (in quotations) to the mapMog folder to access utility data required to execute the function.")}
  if(substr(source.path,nchar(source.path),nchar(source.path)) %in% "/"){source.path <- substr(source.path,0,nchar(source.path)-1)} #remove trailing slash from source path, if present
  if(length(list.dirs(source.path))==0){stop(paste("source.path is not valid. Double check that the file path exists, and that you use forward slashes (/) instead of back slashes (\\)."))}
  if(interpolate == TRUE & is.na(resolution)){stop("a spatial resolution must be provided to interpolate (interploate = TRUE). The suggested resolution is 1590 m (~1 mile).")}
  if(interpolate == TRUE & resolution == 0){stop("spatial resolution must be greater than 1. The suggested resolution is 1590 m (~1 mile).")}
  if(!(round(resolution/30) == resolution/30)){stop("resolution must be a whole number divisible by 30 due to the resolution of underlying datasets. The suggested resolution is 1590 meters (~1 mile).")}
  if(!(round(resolution) == resolution)){stop("resolution must be a whole number divisible by 30 due to the resolution of underlying datasets. The suggested resolution is 1590 meters (~1 mile).")}
  if(interpolate == TRUE & !(power %in% c(1,2,3))){stop("power value must be either 1, 2, or 3 for interpolation (interpolate = TRUE).")}
  if(interpolate == TRUE){
    
    canopy.path <- list.files(paste(source.path,"/utility_canopy",sep=""), #load whatever raster is in the utility_canopy folder
                              pattern = "*.tif$",
                              full.names = TRUE)
    
    if(length(canopy.path)==0){stop("a canopy cover raster is required for interpolation - download NLCD canopy layer (tif) and add it to the 'utility_canopy' in your sourcepath folder.")}
    
  }
  #load/install necessary packages
  if (!require(sf)) message("Installing required package 'sf.'")
  if (!require(sf)) install.packages('sf') #Does user have the package downloaded? If so move on, if not, download.
  library(sf) #load library
  
  if (!require(terra)) message("Installing required package 'terra.'")
  if (!require(terra)) install.packages('terra') #Does user have the package downloaded? If so move on, if not, download.
  library(terra) #load library
  
  #reproject locale to albers equal area conic
  locale <- st_transform(x = locale,
                         crs = 9822) #albers
  
  #process data using rFIA package
  if(API %in% c("rfia","rFIA","RFIA","Rfia","RFia")){ #begin instructions for if rFIA package is being used (allowing for different spellings)
    
    #install/load packages specific to this approach
    if (!require(rFIA)) message("Installing required package 'rFIA.'")
    if (!require(rFIA)) devtools::install_github('hunter-stanke/rFIA') #Does user have the package downloaded? If so move on, if not, download.
    library(rFIA) #load library
    
    
    #determine the location of locale
    message("Extracting locale information") #print progress message
    USA <- st_read(paste(source.path,"/utility_USA/tl_2020_us_state.shp",sep=""))
    USA <- st_transform(x = USA,
                        crs = st_crs(locale)) 
    
    locale.state <- st_filter(x = USA, #clip state to the locale to identify the state in which the locale resides
                              y = st_centroid(x = st_geometry(locale),
                                              of_largest_polygon = TRUE)) #implement a negative buffer for clipping to avoid locales that share jurisdictional borders accidentally crossing due to low-resolution geographic data
    
    #throw errors for locations that are not supported
    if(nrow(locale.state) > 1){stop("Multi-state processing is not supported - run this function seperately for each state-specific component")}
    if(nrow(locale.state) == 0){stop("Functionality for locales outside of the Contiguous United States is not supported")}
    if(locale.state$STUSPS %in% c("DC")){stop("Standard FIA surveys are not performed in the District of Columbia.")}
    
    #extract FIA plots
    message("Extracting FIA plot data at the state level") #print progress message
    options(timeout=3600) #allow downloads to take up to an hour (default in R is 60 seconds)
    if(!(locale.state$STUSPS %in% c("CA","WY","SD","NE","KS","CO","UT","NV","ID","WA","OR","MT","SD"))){ #begin FIA data fetch for all regions which only need the standard tables
      FIA.data <- getFIA(states = locale.state$STUSPS, #use state names isolated above to specify geographic range of FIA data
                         tables = c("PLOT","COND","TREE")) #pull plot, condition, tree, and data only 
    } #end FIA data fetch for all regions which only need the standard tables
    
    if(locale.state$STUSPS %in% c("CA","NV")){#begin FIA data fetch for Region 5, which needs the sitetree table
      FIA.data <- getFIA(states = locale.state$STUSPS, #use state names isolated above to specify geographic range of FIA data
                         tables = c("PLOT","COND","TREE","PLOTGEOM","DWM_COARSE_WOODY_DEBRIS")) #pull plot, condition, tree, and sitetree data only 
      
      FIA.geomplot <- FIA.data$PLOTGEOM
      FIA.woodyDebris <- FIA.data$DWM_COARSE_WOODY_DEBRIS
      
      message("Loading landscape data specific to the Pacific regions")
      paz <- rast(paste(source.path,"/utility_R5_PAZ.tif", sep="")) #load plant association zone (spat)raster
      nwfp <- st_read(paste(source.path,"/utility_NWFPboundary/NWFPboundary.shp", sep="")) #load boundary for northwest forest plan
      nwfp <- st_transform(x = nwfp,
                           crs = st_crs(locale))
    }#endFIA data fetch for Region 5
    
    if(locale.state$STUSPS %in% c("WY","SD","NE","KS","CO","UT","ID","MT","ND")){#begin FIA data fetch for Regions 1, 2 and 4, which may need the plotgeom and vegetation tables
      FIA.data <- getFIA(states = locale.state$STUSPS, #use state names isolated above to specify geographic range of FIA data
                         tables = c("PLOT","COND","TREE","PLOTGEOM","P2VEG_SUBPLOT_SPP")) #pull plot, condition, tree, vegetation, and plotgeom data only 
      FIA.geomplot <- FIA.data$PLOTGEOM
      FIA.veg <- FIA.data$P2VEG_SUBPLOT_SPP
      
      
    }#endFIA data fetch for Regions 1, 2 and 4
    
    if(locale.state$STUSPS %in% c("MT","ID","ND","WY","SD","WA","OR")){
      
      message("Loading landscape data specific to the Northern and Pacific regions")
      ContDivideEast <- st_read(paste(source.path,"/utility_MTcontDivide/utility_MTcontDivide.shp", sep=""))
      ContDivideEast <- st_transform(x = ContDivideEast,
                                     crs = st_crs(locale))
      
      USFStrees <- read.csv(paste(source.path,"/utility.v9-5_2024-10_Natl_MasterTreeSpeciesList.csv",sep=""))
      
      paz <- rast(paste(source.path,"/utility_R5_PAZ.tif", sep="")) #load plant association zone (spat)raster
      
      nwfp <- st_read(paste(source.path,"/utility_NWFPboundary/NWFPboundary.shp", sep="")) #load boundary for northwest forest plan
      nwfp <- st_transform(x = nwfp,
                           crs = st_crs(locale))
      
      counties <- st_read(paste(source.path,"/utility_ORcounties/counties.shp", sep=""))
      counties <- st_transform(x = counties,
                               crs = st_crs(locale))
      
    }
    
    if(locale.state$STUSPS %in% c("WA","OR")){#begin FIA data fetch for Region 6, which needs the plant association zone raster
      FIA.data <- getFIA(states = locale.state$STUSPS, #use state names isolated above to specify geographic range of FIA data
                         tables = c("PLOT","COND","TREE","PLOTGEOM","DWM_COARSE_WOODY_DEBRIS","P2VEG_SUBPLOT_SPP")) #pull plot, condition, tree, and sitetree data only 
      FIA.geomplot <- FIA.data$PLOTGEOM
      FIA.woodyDebris <- FIA.data$DWM_COARSE_WOODY_DEBRIS
      FIA.veg <- FIA.data$P2VEG_SUBPLOT_SPP
      
    }#endFIA data fetch for Region 6
    
    
    #isolate, thin, and convert plot data to sf format
    FIA.plots <- FIA.data$PLOT
    if(recent == "TRUE"){ FIA.plots <- FIA.plots[FIA.plots$INVYR >= as.numeric(substr(as.character(Sys.Date()),1,4)) - 10, ]} #thin to only plots measured in the last ten years, if recent is selected
    FIA.plots <- st_as_sf(x = FIA.plots, 
                          coords = c("LON","LAT"),
                          crs = 4326) #WGS84 (original crs)
    FIA.plots <- st_transform(x = FIA.plots,
                              crs = 9822) #reproject to Albers Equal Area
    
    #clip FIA plots to within 1,000m of the locale to capture FIA plots just outside the boundary that may be helpful in imputation
    FIA.plots <- st_filter(x = FIA.plots,
                           y = st_union(st_buffer(x = locale,
                                                  dist = 5000)))
    
    #isolate condition and tree data corresponding to thinned plots
    FIA.condition <- FIA.data$COND
    FIA.trees <- FIA.data$TREE
    
    time.vector <- vector() #create an empty vector to fill with processing times 
    northern.trigger <- vector() #create an empty vector to fill if northern data are processed 
    
    #assess conditions at the plot level
    for(p in 1:nrow(FIA.plots)){#begin looping through plots (p)
      
      if(p <= 10){start.time <- Sys.time()}
      
      local.plot <- FIA.plots[p, ] #isolate local plot
      MOG.vector <- c(0) #create an vector to fill with MOG statuses since some forests have multiple definitions (but assign the default value of zero to indicate non-MOG)
      
      #trim condition and tree data to the plot level
      conditions <- FIA.condition[which(FIA.condition$PLOT %in% local.plot[1,"PLOT"]), ]
      conditions <- conditions[which(conditions$COND_STATUS_CD == 1), ] #thin to only condition class 1, which is accessible forested land
      conditions <- conditions[which(conditions$INVYR >= substr(Sys.Date()-(365*15),1,4)), ] #retain only plots that have been measured in the last 15 years to avoid pulling data from retired sites
      
      if(nrow(conditions)>0){ #begin instructions if conditions are present
        most.recent.survey <- max(conditions$INVYR) #identify the most recent year the plot was surveyed
        conditions <- conditions[which(conditions$INVYR %in% most.recent.survey), ] #keep only records from the most recent survey
        conditions <- conditions[which(!(is.na(conditions$FLDTYPCD))), ] #remove any condition for which forest type is not specified
        
        condition.issue <- vector() #create an empty vector to fill with data indicating an issue with the condition (saves erroneous processing issues downstream)
        
        for(c in 1:nrow(conditions)){#begin looping through conditions (c) 
          
          #isolate to data of interest
          local.condition <- conditions[c, ] #isolate condition 
          local.trees <- FIA.trees[which(FIA.trees$PLT_CN %in% local.plot[1,"CN"]), ] #retain only trees within the plot
          local.trees <- local.trees[which(local.trees$CONDID %in% local.condition[1,"CONDID"]), ] #retain only trees within the condition (within the plot)
          
          if(nrow(local.trees) == 0){condition.issue <- c(condition.issue, 1)} #add 1 to the condition issue vector to catch errors
          
          if(nrow(local.trees)>0){#begin instructions for if there are trees in the condition
            
            #determine plot region
            region <- ifelse(substr(local.condition$ADFORCD,1,1)[1] %in% "1", "northern",
                             ifelse(substr(local.condition$ADFORCD,1,1)[1] %in% "2", "rocky",
                                    ifelse(substr(local.condition$ADFORCD,1,1)[1] %in% "3", "southwest",
                                           ifelse(substr(local.condition$ADFORCD,1,1)[1] %in% "4", "intermountain",
                                                  ifelse(substr(local.condition$ADFORCD,1,1)[1] %in% "5", "pacific southwest",
                                                         ifelse(substr(local.condition$ADFORCD,1,1)[1] %in% "6", "pacific northwest",
                                                                ifelse(substr(local.condition$ADFORCD,1,1)[1] %in% "8", "southern",
                                                                       ifelse(substr(local.condition$ADFORCD,1,1)[1] %in% "9", "eastern",NA))))))))
            
            
            if(is.na(region)){#begin instructions if plot region is not known
              
              #determine which state the FIA plot is in
              point.intersect <- st_filter(x = USA,
                                           y = local.plot)
              local.state <- point.intersect$STUSPS 
              
              #determine plot region (this isn't perfect, but is accurate for most states)
              region <- ifelse(local.state %in% "MT", "northern",
                               ifelse(local.state %in% c("CO","KA","NE"), "rocky",
                                      ifelse(local.state %in% c("AZ","NM"), "southwest",
                                             ifelse(local.state %in% c("NV","UT"), "intermountain",
                                                    ifelse(local.state %in% c("TX","OK","AR","LA","MS","KY","TN","AL","GA","FL","SC","NC","VA"), "southern",
                                                           ifelse(local.state %in% "CA", "pacific southwest",
                                                                  ifelse(local.state %in% c("WA","OR"), "pacific northwest",
                                                                         ifelse(local.state %in% c("MN","IA","MO","IL","WI","IN","MI","OH","WV","DC","PA","MD","DE","NJ","CT","RI","CT","MA","ME","NH","VT","NY"), "eastern",NA))))))))
              
            }#end instructions if plot region is not known
            
            if(region %in% "northern"){northern.trigger <- c(northern.trigger,1)} #note that northern data was processed to trigger a specific warning message
            
            #determine forest type
            forest.type <- local.condition$FLDTYPCD #see FIADB User Guide Appendix D for forest group codes
            #each of the four points spread around the central plot point has a radius of 24 feet, resulting in an area of 0.0417 acres per condition. Large tree sampling occurs within these 24 foot plots.
            #the unadjusted proportion (CONDPROP_UNADJ) describes the proportion of the total plot comprised by this classification. Use this to back-calculate how many subplots are included in the classification.
            #see section 1.2.5 of FIADB User Guide
            sub.plot.area <- ifelse(local.condition$PROP_BASIS %in% "MACR", 0.2460, 0.0417) #if the sampling took place in a macroplot, use area of macroplot. Otherwise, use regular subplot area (in acres)
            condition.area <- sub.plot.area * (4*local.condition$CONDPROP_UNADJ) #calculate the condition area in acres
            
            #determine stand age
            stand.age <- max(local.condition$STDAGE, local.condition$FLDAGE, na.rm = TRUE) #determine the oldest tree in the stand
            if(stand.age == 999){stand.age <- local.condition$STDAGE} #replace numeric value for no data (999) if that is a problem
            raw.stand.age <- stand.age #save raw stand age for later use (stand age will be turned to binary and will become unrecognizable)
            
            #create mature dataframe for several indices
            mat.df <- local.trees[which(local.trees$DIA >= 1), ] #retain only trees with a diameter greater than 1 inch
            mat.df <- mat.df[which(mat.df$STATUSCD == 1), ] #retain only living trees
            mat.df <- mat.df[which(mat.df$CCLCD %in% c(1,2,3)), ] #retain only dominant trees
            
            if(nrow(mat.df) > 0){#begin calculating indices/metrics if data are available
              
              #calculate variables needed to assess mature status, rather than old growth (based on Table 3 of MOG document)
              tpadom <- nrow(mat.df)/condition.area #calculate dominant live trees per acre
              
              badom <- sum(mat.df$TPA_UNADJ*pi*(mat.df$DIA / 24)*2) #calculate total basal area of dominant trees
              
              QMDdom <- sqrt(badom / (tpadom * 0.005454)) #calculate quadratic mean diameter of dominant trees
              
              local.trees$TPA_class <- ifelse(local.trees$DIA >= 2 & local.trees$DIA <= 9.8, "class 0", #classify trees into diameter classes for Diameter Diversity Index (USFS link to explain this calculation was broken, so this approach is from artificial intelligence that I'm hoping provides a correct summary of the concept)
                                              ifelse(local.trees$DIA >= 9.9 & local.trees$DIA <= 19.7, "class 1",
                                                     ifelse(local.trees$DIA >= 19.8 & local.trees$DIA <= 39.4, "class 2",
                                                            ifelse(local.trees$DIA >= 39.5, "class 3", NA))))
              ddi.vector <- as.vector(local.trees$TPA_class) #convert diameter classes to a vector
              ddiscore <- ((length(ddi.vector[which(ddi.vector %in% "class 0")])/length(ddi.vector))^2) + ((length(ddi.vector[which(ddi.vector %in% "class 1")])/length(ddi.vector))^2) + ((length(ddi.vector[which(ddi.vector %in% "class 2")])/length(ddi.vector))^2) + ((length(ddi.vector[which(ddi.vector %in% "class 3")])/length(ddi.vector))^2) #calculate sum of squared proportion of diameter classes
              ddiscore <- 1 - ddiscore #substract Diameter Diversity from 1 to complete the calculation  
              
              height.vector <- vector() #create an empty vector to fill
              for(k in 1:nrow(local.trees)){#begin looping through local trees (k)
                if(!is.na(local.trees[k,"STATUSCD"]) & local.trees[k,"STATUSCD"] == 1){#begin instructions if the tree is living
                  local.vector <- rep(x = local.trees[k,"HT"], #repeat tree's height for each tree that this point represents in the plot
                                      times = floor(local.trees[k,"TPA_UNADJ"]))
                  
                  height.vector <- c(height.vector, local.vector) #add weighted height record to larger list of tree heights
                } #end instructions if the tree is living
              }#end looping through local trees (k)
              
              top.25prcnt <- quantile(x = height.vector, #identify the cutoff for the top 25% tallest trees in the sample
                                      probs = 0.75,
                                      na.rm=TRUE)
              HTquart <- mean(height.vector[height.vector >= top.25prcnt]) #calculate the mean height of the tallest 25% of trees
              
              HTsd <- sd(height.vector) #calculate the standard deviation of height for all (TPA-weighted) trees
              
              local.dead.trees <- local.trees[which(local.trees$STATUSCD == 2), ] #retain only dead trees
              snagbatot <- sum(local.dead.trees$TPA_UNADJ*pi*(local.dead.trees$DIA / 24)*2, na.rm = TRUE) #calculate total basal area of dead trees
              
              local.MOG.status <- 0 #assign MOG status as 0 until it can be proven to be mature or old growth
              
              if(sum(is.na(c(snagbatot,HTsd,HTquart,ddiscore,QMDdom,badom,tpadom))) == 0){#begin processing if all tree metrics are available
                if(region %in% "eastern"){#begin instructions for eastern region 
                  
                  #assess mature and old growth status based on Tables 17 and 19 in MOG document
                  if(forest.type %in% c(805)){#begin instructions for (eastern) beech maple basswood
                    
                    #determine if information meets MOG criteria
                    stand.age <- ifelse(stand.age >= 140, 1, 0)
                    
                    large.trees <- local.trees[which(local.trees$DIA >= 16), ] #isolate large trees per the forest type definition of large
                    large.trees <- large.trees[which(large.trees$STATUSCD == 1), ] #isolate to living trees only 
                    if(nrow(large.trees) > 0){live.large.tree.dens <- nrow(large.trees)/condition.area}
                    if(nrow(large.trees) == 0){live.large.tree.dens <- 0}
                    live.large.tree.dens <- ifelse(live.large.tree.dens >= 10, 1, 0)
                    
                    local.MOG.status <- ifelse(sum(stand.age,live.large.tree.dens,na.rm=TRUE) == 2, 1, 0) #if both conditions are met, assign 1 for MOG. Otherwise assign 0 for non-MOG
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                  }#end instructions for (eastern) beech maple basswood
                  
                  if(forest.type %in% c(520, 801, 802, 809)){#begin instructions for (eastern) northern hardwood
                    
                    #determine if information meets MOG criteria
                    stand.age <- ifelse(stand.age >= 140, 1, 0)
                    
                    large.trees <- local.trees[which(local.trees$DIA >= 16), ] #isolate large trees per the forest type definition of large
                    large.trees <- large.trees[which(large.trees$STATUSCD == 1), ] #isolate to living trees only 
                    if(nrow(large.trees) > 0){live.large.tree.dens <- nrow(large.trees)/condition.area}
                    if(nrow(large.trees) == 0){live.large.tree.dens <- 0}
                    live.large.tree.dens <- ifelse(live.large.tree.dens >= 10, 1, 0)
                    
                    local.MOG.status <- ifelse(sum(stand.age,live.large.tree.dens,na.rm=TRUE) == 2, 1, 0) #if both conditions are met, assign 1 for MOG. Otherwise assign 0 for non-MOG
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end instructions for (eastern) northern hardwood
                  
                  if(forest.type %in% c(162, 163, 165, 167, 182, 184, 404, 405, 501, 502, 506, 507, 509, 510, 513, 515)){#begin instructions for (eastern) dry oak
                    
                    #determine if information meets MOG criteria
                    stand.age <- ifelse(stand.age >= 100, 1, 0)
                    
                    large.trees <- local.trees[which(local.trees$DIA >= 16), ] #isolate large trees per the forest type definition of large
                    large.trees <- large.trees[which(large.trees$STATUSCD == 1), ] #isolate to living trees only 
                    if(nrow(large.trees) > 0){live.large.tree.dens <- nrow(large.trees)/condition.area}
                    if(nrow(large.trees) == 0){live.large.tree.dens <- 0}
                    live.large.tree.dens <- ifelse(live.large.tree.dens >= 20, 1, 0)
                    
                    local.MOG.status <- ifelse(sum(stand.age,live.large.tree.densna.rm=TRUE) == 2, 1, 0) #if both conditions are met, assign 1 for MOG. Otherwise assign 0 for non-MOG
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end instructions for (eastern) dry oak
                  
                  if(forest.type %in% c(503, 504, 505, 511, 512, 516)){#begin instructions for (eastern) mesic northern oak
                    
                    #determine if information meets MOG criteria
                    stand.age <- ifelse(stand.age >= 160, 1, 0)
                    
                    large.trees <- local.trees[which(local.trees$DIA >= 20), ] #isolate large trees per the forest type definition of large
                    large.trees <- large.trees[which(large.trees$STATUSCD == 1), ] #isolate to living trees only 
                    if(nrow(large.trees) > 0){live.large.tree.dens <- nrow(large.trees)/condition.area}
                    if(nrow(large.trees) == 0){live.large.tree.dens <- 0}
                    live.large.tree.dens <- ifelse(live.large.tree.dens >= 5, 1, 0)
                    
                    local.MOG.status <- ifelse(sum(stand.age,live.large.tree.dens,na.rm=TRUE) == 2, 1, 0) #if both conditions are met, assign 1 for MOG. Otherwise assign 0 for non-MOG
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end instructions for (eastern) mesic northern oak
                  
                  if(forest.type %in% c(701, 702, 703, 704, 705, 706, 707, 708, 709)){#begin instructions for (eastern) wetland hardwood
                    
                    #determine if information meets MOG criteria
                    stand.age <- ifelse(stand.age >= 120, 1, 0)
                    
                    large.trees <- local.trees[which(local.trees$DIA >= 18), ] #isolate large trees per the forest type definition of large
                    large.trees <- large.trees[which(large.trees$STATUSCD == 1), ] #isolate to living trees only 
                    if(nrow(large.trees) > 0){live.large.tree.dens <- nrow(large.trees)/condition.area}
                    if(nrow(large.trees) == 0){live.large.tree.dens <- 0}
                    live.large.tree.dens <- ifelse(live.large.tree.dens >= 10, 1, 0)
                    
                    local.MOG.status <- ifelse(sum(stand.age,live.large.tree.dens,na.rm=TRUE) == 2, 1, 0) #if both conditions are met, assign 1 for MOG. Otherwise assign 0 for non-MOG
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end instructions for (eastern) wetland hardwood
                  
                  if(forest.type %in% c(104, 105, 401)){#begin instructions for (eastern) conifer northern hardwood
                    
                    #determine if information meets MOG criteria
                    stand.age <- ifelse(stand.age >= 140, 1, 0)
                    
                    large.trees <- local.trees[which(local.trees$DIA >= 16), ] #isolate large trees per the forest type definition of large
                    large.trees <- large.trees[which(large.trees$STATUSCD == 1), ] #isolate to living trees only 
                    if(nrow(large.trees) > 0){live.large.tree.dens <- nrow(large.trees)/condition.area}
                    if(nrow(large.trees) == 0){live.large.tree.dens <- 0}
                    live.large.tree.dens <- ifelse(live.large.tree.dens >= 10, 1, 0)
                    
                    local.MOG.status <- ifelse(sum(stand.age,live.large.tree.dens,na.rm=TRUE) == 2, 1, 0) #if both conditions are met, assign 1 for MOG. Otherwise assign 0 for non-MOG
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end instructions for (eastern) conifer northern hardwood
                  
                  if(forest.type %in% c(101, 102, 103)){#begin instructions for (eastern) northern pine
                    
                    #determine if information meets MOG criteria
                    stand.age <- ifelse(stand.age >= 100, 1, 0)
                    
                    large.trees <- local.trees[which(local.trees$DIA >= 12), ] #isolate large trees per the forest type definition of large
                    large.trees <- large.trees[which(large.trees$STATUSCD == 1), ] #isolate to living trees only 
                    if(nrow(large.trees) > 0){live.large.tree.dens <- nrow(large.trees)/condition.area}
                    if(nrow(large.trees) == 0){live.large.tree.dens <- 0}
                    live.large.tree.dens <- ifelse(live.large.tree.dens >= 20, 1, 0)
                    
                    local.MOG.status <- ifelse(sum(stand.age,live.large.tree.dens,na.rm=TRUE) == 2, 1, 0) #if both conditions are met, assign 1 for MOG. Otherwise assign 0 for non-MOG
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end instructions for (eastern) northern pine
                  
                  if(forest.type %in% c(121, 123, 124, 128, 129)){#begin instructions for (eastern) montane spruce
                    
                    #determine if information meets MOG criteria
                    stand.age <- ifelse(stand.age >= 140, 1, 0)
                    
                    large.trees <- local.trees[which(local.trees$DIA >= 15), ] #isolate large trees per the forest type definition of large
                    large.trees <- large.trees[which(large.trees$STATUSCD == 1), ] #isolate to living trees only 
                    if(nrow(large.trees) > 0){live.large.tree.dens <- nrow(large.trees)/condition.area}
                    if(nrow(large.trees) == 0){live.large.tree.dens <- 0}
                    live.large.tree.dens <- ifelse(live.large.tree.dens >= 10, 1, 0)
                    
                    local.MOG.status <- ifelse(sum(stand.age,live.large.tree.dens,na.rm=TRUE) == 2, 1, 0) #if both conditions are met, assign 1 for MOG. Otherwise assign 0 for non-MOG
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end instructions for (eastern) montane spruce
                  
                  if(forest.type %in% c(122, 125)){#begin instructions for (eastern) Sub-boreal spruce/fir
                    
                    #determine if information meets MOG criteria
                    stand.age <- ifelse(stand.age >= 140, 1, 0)
                    
                    large.trees <- local.trees[which(local.trees$DIA >= 12), ] #isolate large trees per the forest type definition of large
                    large.trees <- large.trees[which(large.trees$STATUSCD == 1), ] #isolate to living trees only 
                    if(nrow(large.trees) > 0){live.large.tree.dens <- nrow(large.trees)/condition.area}
                    if(nrow(large.trees) == 0){live.large.tree.dens <- 0}
                    live.large.tree.dens <- ifelse(live.large.tree.dens >= 10, 1, 0)
                    
                    local.MOG.status <- ifelse(sum(stand.age,live.large.tree.dens,na.rm=TRUE) == 2, 1, 0) #if both conditions are met, assign 1 for MOG. Otherwise assign 0 for non-MOG
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end instructions for (eastern) Sub-boreal spruce/fir
                  
                  if(!(forest.type %in% c(122, 125, 121, 123, 124, 128, 129, 101, 102, 103, 104, 105, 401, 701, 702, 703, 704, 705, 706, 707, 708, 709, 503, 504, 505, 511, 512, 516, 162, 163, 165, 167, 182, 184, 404, 405, 501, 502, 506, 507, 509, 510, 513, 515, 520, 801, 802, 809, 805))){#begin instructions for (eastern) all other forest types
                    
                    #determine if information meets MOG criteria
                    stand.age <- ifelse(stand.age >= 100, 1, 0)
                    
                    large.trees <- local.trees[which(local.trees$DIA >= 14), ] #isolate large trees per the forest type definition of large
                    large.trees <- large.trees[which(large.trees$STATUSCD == 1), ] #isolate to living trees only 
                    if(nrow(large.trees) > 0){live.large.tree.dens <- nrow(large.trees)/condition.area}
                    if(nrow(large.trees) == 0){live.large.tree.dens <- 0}
                    live.large.tree.dens <- ifelse(live.large.tree.dens >= 10, 1, 0)
                    
                    local.MOG.status <- ifelse(sum(stand.age,live.large.tree.dens,na.rm=TRUE) == 2, 1, 0) #if both conditions are met, assign 1 for MOG. Otherwise assign 0 for non-MOG
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end instructions for (eastern) all other forest types
                  
                  #if not old growth, assess to see if the plot is at least mature
                  if(local.MOG.status == 0){#begin mature assessment
                    
                    if(forest.type %in% c(104,105,401,400:409)){#begin mature assessment for (eastern) conifer northern hardwood
                      
                      t.QMDdom <- ifelse(QMDdom >= 14, 0.3, 0)
                      t.badom <- ifelse(badom >= 104.3, 0.22, 0)
                      t.snagbatot <- ifelse(snagbatot >= 14.5, 0.19, 0)
                      t.tpadom <- ifelse(tpadom <= 73.4, 0.16, 0)
                      t.HTsd <- ifelse(HTsd >= 32, 0.13, 0)
                      
                      local.MOG.status <- sum(t.QMDdom,t.badom,t.snagbatot,t.tpadom,t.HTsd,na.rm=TRUE) #sum maturity thresholds to achieve maturity index (0.5 is the cutoff)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for (eastern) conifer northern hardwood
                    
                    if(forest.type %in% c(520,801,802,809,900:905,805,600:609,500:520,960:962,701,702,703,704,705,706,707,708,709)){#begin mature assessment for (eastern) northern hardwood
                      
                      t.QMDdom <- ifelse(QMDdom >= 9.9, 0.29, 0)
                      t.HTquart <- ifelse(HTquart >= 43.3, 0.2, 0)
                      t.badom <- ifelse(badom >= 60.9, 0.18,0)
                      t.tpadom <- ifelse(tpadom <= 97.6, 0.18, 0)
                      t.HTsd <- ifelse(HTsd >= 32.9, 0.14, 0)
                      
                      local.MOG.status <- sum(t.QMDdom,t.HTquart,t.badom,t.tpadom,t.HTsd,na.rm=TRUE)#sum maturity thresholds to achieve maturity index (0.5 is the cutoff)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for (eastern) northern hardwood
                    
                    if(forest.type %in% c(101,102,103,160:168,380:385,390,391)){#begin mature assessment for (eastern) northern pine
                      
                      t.QMDdom <- ifelse(QMDdom >= 11.9, 0.3, 0)
                      t.HTsd <- ifelse(HTsd >= 67.4, 0.22, 0)
                      t.HTquart <- ifelse(HTquart >= 38, 0.21, 0)
                      t.tpadom <- ifelse(tpadom <= 83.2, 0.18, 0)
                      t.badom <- ifelse(badom >- 81.5, 0.09, 0)
                      
                      local.MOG.status <- sum(t.QMDdom,t.HTsd,t.HTquart,t.tpadom,t.badom,na.rm=TRUE) #sum maturity thresholds to achieve maturity index (0.5 is the cutoff)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for (eastern) northern pine
                    
                    if(forest.type %in% c(162,173,165,167,182,184,404,405,501,502,506,507,509,510,513,515,503,504,505,511,512,516)){#begin mature assessment for (eastern) oak
                      
                      t.QMDdom <- ifelse(QMDdom >= 12.7, 0.37, 0)
                      t.tpadom <- ifelse(tpadom <= 73.4, 0.22, 0)
                      t.HTquart <- ifelse(HTquart >= 52.9, 0.21, 0)
                      t.HTsd <- ifelse(HTsd >= 36.5, 0.2, 0)
                      
                      local.MOG.status <- sum(t.QMDdom,t.tpadom,t.HTquart,t.HTsd,na.rm=TRUE) #sum maturity thresholds to achieve maturity index (0.5 is the cutoff)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for (eastern) oak
                    
                    if(forest.type %in% c(122,125,121,123,124,128,129,120:129)){#begin mature assessment for spruce/fir
                      
                      t.ddiscore <- ifelse(ddiscore >= 22.2, 0.4, 0)
                      t.badom <- ifelse(badom >= 76.2, 0.36, 0)
                      t.HTquart <- ifelse(HTquart >= 32, 0.24, 0)
                      
                      local.MOG.status <- sum(t.ddiscore,t.badom,t.HTquart,na.rm=TRUE)#sum maturity thresholds to achieve maturity index (0.5 is the cutoff)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for spruce/fir
                    
                  }#end mature assessment
                  
                }#end instructions for eastern region
                
                if(region %in% "southern"){#begin instructions for southern region
                  
                  #identify applicable southern forest types (crosswalk) since southern region can fall into multiple forest categories (based on Table 16 of MOG document)
                  southern.type <- vector() #create an empty vector to fill with forest types
                  if(forest.type %in% c(104,105,123,124)){southern.type <- c(southern.type,2)}
                  if(forest.type %in% c(129)){southern.type <- c(southern.type,31)}
                  if(forest.type %in% c(141)){southern.type <- c(southern.type,26,29)}
                  if(forest.type %in% c(142,166,407)){southern.type <- c(southern.type,29)}
                  if(forest.type %in% c(161)){southern.type <- c(southern.type,25)}
                  if(forest.type %in% c(162,163,404,405,409)){southern.type <- c(southern.type,24,25)}
                  if(forest.type %in% c(165,167)){southern.type <- c(southern.type,24)}
                  if(forest.type %in% c(400)){southern.type <- c(southern.type,2,24,25,26,29)}
                  if(forest.type %in% c(401)){southern.type <- c(southern.type,2)}
                  if(forest.type %in% c(403)){southern.type <- c(southern.type,26)}
                  if(forest.type %in% c(406)){southern.type <- c(southern.type,25)}
                  if(forest.type %in% c(500)){southern.type <- c(southern.type,5,13,21,22,24,27)}
                  if(forest.type %in% c(501)){southern.type <- c(southern.type,22)}
                  if(forest.type %in% c(502,515,519)){southern.type <- c(southern.type,21,22)}
                  if(forest.type %in% c(504)){southern.type <- c(southern.type,21,27)}
                  if(forest.type %in% c(505)){southern.type <- c(southern.type,21)}
                  if(forest.type %in% c(506,511,516)){southern.type <- c(southern.type,5)}
                  if(forest.type %in% c(508)){southern.type <- c(southern.type,13)}
                  if(forest.type %in% c(510)){southern.type <- c(southern.type,21,22,24)}
                  if(forest.type %in% c(514)){southern.type <- c(southern.type,22,24)}
                  if(forest.type %in% c(517,800,801,805)){southern.type <- c(southern.type,1,5)}
                  if(forest.type %in% c(520)){southern.type <- c(southern.type,27)}
                  if(forest.type %in% c(600)){southern.type <- c(southern.type,6,10,13,22,27,28)}
                  if(forest.type %in% c(601,602,605,706)){southern.type <- c(southern.type,13)}
                  if(forest.type %in% c(607,609)){southern.type <- c(southern.type,14)}
                  if(forest.type %in% c(608,809)){southern.type <- c(southern.type,10)}
                  if(forest.type %in% c(700)){southern.type <- c(southern.type,10,28)}
                  if(forest.type %in% c(702,703,704)){southern.type <- c(southern.type,28)}
                  if(forest.type %in% c(705)){southern.type <- c(southern.type,13,28)}
                  if(forest.type %in% c(708)){southern.type <- c(southern.type,10,13)}
                  if(forest.type %in% c(709)){southern.type <- c(southern.type,28)}
                  if(forest.type %in% c(902)){southern.type <- c(southern.type,31,2)} 
                  if(forest.type %in% c(962)){southern.type <- c(southern.type,1,5,6,10,13,21,22,27,28)}
                  
                  if(length(southern.type > 0)){#begin instructions if condition fits within a defined southern forest type
                    #southern.MOG.status <- vector() #create an empty vector to fill with MOG status based on each forest type the FIA plot falls under
                    for(s in 1:length(southern.type)){#begin looping through southern forest types (s)
                      
                      #pull old growth qualifying information for condition (based on table 15 of MOG document)
                      stand.age <- max(local.condition$STDAGE, local.condition$FLDAGE, na.rm = TRUE) #determine the oldest tree in the stand
                      raw.stand.age <- stand.age #save the raw stand age for later use
                      dead.trees <- local.trees[which(local.trees$STATUSCD == 2), ] #isolate to dead trees only 
                      dead.tree.dens <- nrow(dead.trees)/condition.area
                      live.trees <- local.trees[which(local.trees$STATUSCD == 1), ] #isolate to live trees only 
                      
                      #calculate stand basal area since FIA plot data performs calculation on all trees greater than 1 in DBH, but regional MOG definition only wants trees greater than 5 DBH
                      live.basal.trees <- live.trees[live.trees$DIA >= 5, ]
                      live.basal.trees$basal <- (live.basal.trees$DIA^2)*0.005454 #use foresters constant to convert DBH to square feet
                      live.basal.trees$basal <- live.basal.trees$basal*(1/condition.area) #multiply by expansion factor to arive at sqft per acre
                      stand.basal.area <- sum(live.basal.trees$basal)
                      
                      if(southern.type[s] %in% 1){ #begin instructions for (southern) northern hardwood forest
                        
                        #determine if information meets MOG criteria
                        stand.age <- ifelse(stand.age >= 100, 1, 0)
                        
                        large.trees <- live.trees[which(live.trees$DIA >= 14), ] #isolate large trees per the forest type definition of large
                        if(nrow(large.trees) > 0){live.large.tree.dens <- nrow(large.trees)/condition.area} #calculate density of large trees
                        if(nrow(large.trees) == 0){live.large.tree.dens <- 0} #handle large tree density if no large trees exist
                        live.large.tree.dens <- ifelse(live.large.tree.dens >= 6, 1, 0) #assign MOG status (1) if condition meets forest type definition
                        
                        stand.basal.area <- ifelse(stand.basal.area >= 40, 1, 0) #assign MOG status (1) if condition meets forest type definition
                        
                        dead.tree.dens <- ifelse(dead.tree.dens >= 13, 1, 0) #assign MOG status (1) if condition meets forest type definition
                        
                        local.MOG.status <- ifelse(sum(stand.age,dead.tree.dens,live.large.tree.dens,stand.basal.areana.rm=TRUE) == 4, 1, 0) #if all conditions are met, assign 1 for old growth. Otherwise assign 0 for non-OG
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end instructions for (southern) northern hardwood forest  
                      
                      if(southern.type[s] %in% 2){ #begin instructions for (southern) conifer-northern hardwood forest
                        
                        #determine if information meets MOG criteria
                        stand.age <- ifelse(stand.age >= 140, 1, 0)
                        
                        large.trees <- live.trees[which(live.trees$DIA >= 20), ] #isolate large trees per the forest type definition of large
                        if(nrow(large.trees) > 0){live.large.tree.dens <- nrow(large.trees)/condition.area} #calculate density of large trees
                        if(nrow(large.trees) == 0){live.large.tree.dens <- 0} #handle large tree density if no large trees exist
                        live.large.tree.dens <- ifelse(live.large.tree.dens >= 6, 1, 0) #assign MOG status (1) if condition meets forest type definition
                        
                        stand.basal.area <- ifelse(stand.basal.area >= 40, 1, 0) #assign MOG status (1) if condition meets forest type definition
                        
                        dead.tree.dens <- ifelse(dead.tree.dens >= 6, 1, 0) #assign MOG status (1) if condition meets forest type definition
                        
                        local.MOG.status <- ifelse(sum(stand.age,dead.tree.dens,live.large.tree.dens,stand.basal.area,na.rm=TRUE) == 4, 1, 0) #if all conditions are met, assign 1 for old growth. Otherwise assign 0 for non-OG
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end instructions for (southern) conifer-northern hardwood forest  
                      
                      if(southern.type[s] %in% 5){ #begin instructions for (southern) mixed mesophytic and western mesophytic forest
                        
                        #determine if information meets MOG criteria
                        stand.age <- ifelse(stand.age >= 140, 1, 0)
                        
                        large.trees <- live.trees[which(live.trees$DIA >= 30), ] #isolate large trees per the forest type definition of large
                        if(nrow(large.trees) > 0){live.large.tree.dens <- nrow(large.trees)/condition.area} #calculate density of large trees
                        if(nrow(large.trees) == 0){live.large.tree.dens <- 0} #handle large tree density if no large trees exist
                        live.large.tree.dens <- ifelse(live.large.tree.dens >= 6, 1, 0) #assign MOG status (1) if condition meets forest type definition
                        
                        stand.basal.area <- ifelse(stand.basal.area >= 40, 1, 0) #assign MOG status (1) if condition meets forest type definition
                        
                        dead.tree.dens <- ifelse(dead.tree.dens >= 4, 1, 0) #assign MOG status (1) if condition meets forest type definition
                        
                        local.MOG.status <- ifelse(sum(stand.age,dead.tree.dens,live.large.tree.dens,stand.basal.area,na.rm=TRUE) == 4, 1, 0) #if all conditions are met, assign 1 for old growth. Otherwise assign 0 for non-OG
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end instructions for (southern) mixed mesophytic and western mesophytic forest
                      
                      
                      if(southern.type[s] %in% 6){ #begin instructions for (southern) coastal plain upland mesic hardwood forest
                        
                        #determine if information meets MOG criteria
                        stand.age <- ifelse(stand.age >= 120, 1, 0)
                        
                        large.trees <- live.trees[which(live.trees$DIA >= 24), ] #isolate large trees per the forest type definition of large
                        if(nrow(large.trees) > 0){live.large.tree.dens <- nrow(large.trees)/condition.area} #calculate density of large trees
                        if(nrow(large.trees) == 0){live.large.tree.dens <- 0} #handle large tree density if no large trees exist
                        live.large.tree.dens <- ifelse(live.large.tree.dens >= 6, 1, 0) #assign MOG status (1) if condition meets forest type definition
                        
                        stand.basal.area <- ifelse(stand.basal.area >= 40, 1, 0) #assign MOG status (1) if condition meets forest type definition
                        
                        dead.tree.dens <- ifelse(dead.tree.dens >= 4, 1, 0) #assign MOG status (1) if condition meets forest type definition
                        
                        local.MOG.status <- ifelse(sum(stand.age,dead.tree.dens,live.large.tree.dens,stand.basal.area,na.rm=TRUE) == 4, 1, 0) #if all conditions are met, assign 1 for old growth. Otherwise assign 0 for non-OG
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end instructions for (southern) coastal plain upland mesic hardwood forest
                      
                      if(southern.type[s] %in% 10){ #begin instructions for (southern) hardwood wetland forest
                        
                        #determine if information meets MOG criteria
                        stand.age <- ifelse(stand.age >= 120, 1, 0)
                        
                        large.trees <- live.trees[which(live.trees$DIA >= 20), ] #isolate large trees per the forest type definition of large
                        if(nrow(large.trees) > 0){live.large.tree.dens <- nrow(large.trees)/condition.area} #calculate density of large trees
                        if(nrow(large.trees) == 0){live.large.tree.dens <- 0} #handle large tree density if no large trees exist
                        live.large.tree.dens <- ifelse(live.large.tree.dens >= 6, 1, 0) #assign MOG status (1) if condition meets forest type definition
                        
                        stand.basal.area <- ifelse(stand.basal.area >= 40, 1, 0) #assign MOG status (1) if condition meets forest type definition
                        
                        dead.tree.dens <- ifelse(dead.tree.dens >= 0, 1, 0) #assign MOG status (1) if condition meets forest type definition
                        
                        local.MOG.status <- ifelse(sum(stand.age,dead.tree.dens,live.large.tree.dens,stand.basal.area,na.rm=TRUE) == 4, 1, 0) #if all conditions are met, assign 1 for old growth. Otherwise assign 0 for non-OG
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end instructions for (southern) hardwood wetland forest
                      
                      if(southern.type[s] %in% 13){ #begin instructions for (southern) river floodplain hardwood forest
                        
                        #determine if information meets MOG criteria
                        stand.age <- ifelse(stand.age >= 100, 1, 0)
                        
                        large.trees <- live.trees[which(live.trees$DIA >= 16), ] #isolate large trees per the forest type definition of large
                        if(nrow(large.trees) > 0){live.large.tree.dens <- nrow(large.trees)/condition.area} #calculate density of large trees
                        if(nrow(large.trees) == 0){live.large.tree.dens <- 0} #handle large tree density if no large trees exist
                        live.large.tree.dens <- ifelse(live.large.tree.dens >= 6, 1, 0) #assign MOG status (1) if condition meets forest type definition
                        
                        stand.basal.area <- ifelse(stand.basal.area >= 40, 1, 0) #assign MOG status (1) if condition meets forest type definition
                        
                        dead.tree.dens <- ifelse(dead.tree.dens >= 0, 1, 0) #assign MOG status (1) if condition meets forest type definition
                        
                        local.MOG.status <- ifelse(sum(stand.age,dead.tree.dens,live.large.tree.dens,stand.basal.area,na.rm=TRUE) == 4, 1, 0) #if all conditions are met, assign 1 for old growth. Otherwise assign 0 for non-OG
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end instructions for (southern) river floodplain hardwood forest
                      
                      if(southern.type[s] %in% 14){ #begin instructions for (southern) cypress-tupelo swamp forest
                        
                        #determine if information meets MOG criteria
                        stand.age <- ifelse(stand.age >= 120, 1, 0)
                        
                        large.trees <- live.trees[which(live.trees$DIA >= 8), ] #isolate large trees per the forest type definition of large
                        if(nrow(large.trees) > 0){live.large.tree.dens <- nrow(large.trees)/condition.area} #calculate density of large trees
                        if(nrow(large.trees) == 0){live.large.tree.dens <- 0} #handle large tree density if no large trees exist
                        live.large.tree.dens <- ifelse(live.large.tree.dens >= 6, 1, 0) #assign MOG status (1) if condition meets forest type definition
                        
                        stand.basal.area <- ifelse(stand.basal.area >= 40, 1, 0) #assign MOG status (1) if condition meets forest type definition
                        
                        dead.tree.dens <- ifelse(dead.tree.dens >= 3, 1, 0) #assign MOG status (1) if condition meets forest type definition
                        
                        local.MOG.status <- ifelse(sum(stand.age,dead.tree.dens,live.large.tree.dens,stand.basal.area,na.rm=TRUE) == 4, 1, 0) #if all conditions are met, assign 1 for old growth. Otherwise assign 0 for non-OG
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end instructions for (southern) cypress-tupelo swamp forest
                      
                      if(southern.type[s] %in% 21){ #begin instructions for (southern) dry-mesic oak forest
                        
                        #determine if information meets MOG criteria
                        stand.age <- ifelse(stand.age >= 130, 1, 0)
                        
                        large.trees <- live.trees[which(live.trees$DIA >= 20), ] #isolate large trees per the forest type definition of large
                        if(nrow(large.trees) > 0){live.large.tree.dens <- nrow(large.trees)/condition.area} #calculate density of large trees
                        if(nrow(large.trees) == 0){live.large.tree.dens <- 0} #handle large tree density if no large trees exist
                        live.large.tree.dens <- ifelse(live.large.tree.dens >= 6, 1, 0) #assign MOG status (1) if condition meets forest type definition
                        
                        stand.basal.area <- ifelse(stand.basal.area >= 40, 1, 0) #assign MOG status (1) if condition meets forest type definition
                        
                        dead.tree.dens <- ifelse(dead.tree.dens >= 26, 1, 0) #assign MOG status (1) if condition meets forest type definition
                        
                        local.MOG.status <- ifelse(sum(stand.age,dead.tree.dens,live.large.tree.dens,stand.basal.area,na.rm=TRUE) == 4, 1, 0) #if all conditions are met, assign 1 for old growth. Otherwise assign 0 for non-OG
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end instructions for (southern) dry-mesic oak forest
                      
                      if(southern.type[s] %in% 22){ #begin instructions for (southern) dry and xeric oak forest, woodland, and savanna
                        
                        #determine if information meets MOG criteria
                        stand.age <- ifelse(stand.age >= 90, 1, 0)
                        
                        large.trees <- live.trees[which(live.trees$DIA >= 8), ] #isolate large trees per the forest type definition of large
                        if(nrow(large.trees) > 0){live.large.tree.dens <- nrow(large.trees)/condition.area} #calculate density of large trees
                        if(nrow(large.trees) == 0){live.large.tree.dens <- 0} #handle large tree density if no large trees exist
                        live.large.tree.dens <- ifelse(live.large.tree.dens >= 6, 1, 0) #assign MOG status (1) if condition meets forest type definition
                        
                        stand.basal.area <- ifelse(stand.basal.area >= 10, 1, 0) #assign MOG status (1) if condition meets forest type definition
                        
                        dead.tree.dens <- ifelse(dead.tree.dens >= 10, 1, 0) #assign MOG status (1) if condition meets forest type definition
                        
                        local.MOG.status <- ifelse(sum(stand.age,dead.tree.dens,live.large.tree.dens,stand.basal.area,na.rm=TRUE) == 4, 1, 0) #if all conditions are met, assign 1 for old growth. Otherwise assign 0 for non-OG
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end instructions for (southern) dry and xeric oak forest, woodland, and savanna
                      
                      if(southern.type[s] %in% 24){ #begin instructions for (southern) xeric pine and pine-oak forest and woodland
                        
                        #determine if information meets MOG criteria
                        stand.age <- ifelse(stand.age >= 100, 1, 0)
                        
                        large.trees <- live.trees[which(live.trees$DIA >= 10), ] #isolate large trees per the forest type definition of large
                        if(nrow(large.trees) > 0){live.large.tree.dens <- nrow(large.trees)/condition.area} #calculate density of large trees
                        if(nrow(large.trees) == 0){live.large.tree.dens <- 0} #handle large tree density if no large trees exist
                        live.large.tree.dens <- ifelse(live.large.tree.dens >= 6, 1, 0) #assign MOG status (1) if condition meets forest type definition
                        
                        stand.basal.area <- ifelse(stand.basal.area >= 20, 1, 0) #assign MOG status (1) if condition meets forest type definition
                        
                        dead.tree.dens <- ifelse(dead.tree.dens >= 6, 1, 0) #assign MOG status (1) if condition meets forest type definition
                        
                        local.MOG.status <- ifelse(sum(stand.age,dead.tree.dens,live.large.tree.dens,stand.basal.area,na.rm=TRUE) == 4, 1, 0) #if all conditions are met, assign 1 for old growth. Otherwise assign 0 for non-OG
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end instructions for (southern) xeric pine and pine-oak forest and woodland
                      
                      if(southern.type[s] %in% 25){ #begin instructions for (southern) dry and dry-mesic oak-pine forest
                        
                        #determine if information meets MOG criteria
                        stand.age <- ifelse(stand.age >= 120, 1, 0)
                        
                        large.trees <- live.trees[which(live.trees$DIA >= 19), ] #isolate large trees per the forest type definition of large
                        if(nrow(large.trees) > 0){live.large.tree.dens <- nrow(large.trees)/condition.area} #calculate density of large trees
                        if(nrow(large.trees) == 0){live.large.tree.dens <- 0} #handle large tree density if no large trees exist
                        live.large.tree.dens <- ifelse(live.large.tree.dens >= 6, 1, 0) #assign MOG status (1) if condition meets forest type definition
                        
                        stand.basal.area <- ifelse(stand.basal.area >= 40, 1, 0) #assign MOG status (1) if condition meets forest type definition
                        
                        dead.tree.dens <- ifelse(dead.tree.dens >= 15, 1, 0) #assign MOG status (1) if condition meets forest type definition
                        
                        local.MOG.status <- ifelse(sum(stand.age,dead.tree.dens,live.large.tree.dens,stand.basal.area,na.rm=TRUE) == 4, 1, 0) #if all conditions are met, assign 1 for old growth. Otherwise assign 0 for non-OG
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end instructions for (southern) dry and dry-mesic oak-pine forest
                      
                      if(southern.type[s] %in% 26){ #begin instructions for (southern) upland longleaf and south Florida slash pine forest, woodland, and savanna
                        
                        #determine if information meets MOG criteria
                        stand.age <- ifelse(stand.age >= 80, 1, 0)
                        
                        large.trees <- live.trees[which(live.trees$DIA >= 16), ] #isolate large trees per the forest type definition of large
                        if(nrow(large.trees) > 0){live.large.tree.dens <- nrow(large.trees)/condition.area} #calculate density of large trees
                        if(nrow(large.trees) == 0){live.large.tree.dens <- 0} #handle large tree density if no large trees exist
                        live.large.tree.dens <- ifelse(live.large.tree.dens >= 6, 1, 0) #assign MOG status (1) if condition meets forest type definition
                        
                        stand.basal.area <- ifelse(stand.basal.area >= 10, 1, 0) #assign MOG status (1) if condition meets forest type definition
                        
                        dead.tree.dens <- ifelse(dead.tree.dens >= 0, 1, 0) #assign MOG status (1) if condition meets forest type definition
                        
                        local.MOG.status <- ifelse(sum(stand.age,dead.tree.dens,live.large.tree.dens,stand.basal.area,na.rm=TRUE) == 4, 1, 0) #if all conditions are met, assign 1 for old growth. Otherwise assign 0 for non-OG
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end instructions for (southern) upland longleaf and south Florida slash pine forest, woodland, and savanna
                      
                      if(southern.type[s] %in% 27){ #begin instructions for (southern) seasonally wet oak-hardwood woodland
                        
                        #determine if information meets MOG criteria
                        stand.age <- ifelse(stand.age >= 100, 1, 0)
                        
                        large.trees <- live.trees[which(live.trees$DIA >= 20), ] #isolate large trees per the forest type definition of large
                        if(nrow(large.trees) > 0){live.large.tree.dens <- nrow(large.trees)/condition.area} #calculate density of large trees
                        if(nrow(large.trees) == 0){live.large.tree.dens <- 0} #handle large tree density if no large trees exist
                        live.large.tree.dens <- ifelse(live.large.tree.dens >= 6, 1, 0) #assign MOG status (1) if condition meets forest type definition
                        
                        stand.basal.area <- ifelse(stand.basal.area >= 40, 1, 0) #assign MOG status (1) if condition meets forest type definition
                        
                        dead.tree.dens <- ifelse(dead.tree.dens >= 0, 1, 0) #assign MOG status (1) if condition meets forest type definition
                        
                        local.MOG.status <- ifelse(sum(stand.age,dead.tree.dens,live.large.tree.dens,stand.basal.area,na.rm=TRUE) == 4, 1, 0) #if all conditions are met, assign 1 for old growth. Otherwise assign 0 for non-OG
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end instructions for (southern) seasonally wet oak-hardwood woodland
                      
                      if(southern.type[s] %in% 28){ #begin instructions for (southern) eastern riverfront forest
                        
                        #determine if information meets MOG criteria
                        stand.age <- ifelse(stand.age >= 100, 1, 0)
                        
                        large.trees <- live.trees[which(live.trees$DIA >= 25), ] #isolate large trees per the forest type definition of large
                        if(nrow(large.trees) > 0){live.large.tree.dens <- nrow(large.trees)/condition.area} #calculate density of large trees
                        if(nrow(large.trees) == 0){live.large.tree.dens <- 0} #handle large tree density if no large trees exist
                        live.large.tree.dens <- ifelse(live.large.tree.dens >= 6, 1, 0) #assign MOG status (1) if condition meets forest type definition
                        
                        stand.basal.area <- ifelse(stand.basal.area >= 40, 1, 0) #assign MOG status (1) if condition meets forest type definition
                        
                        dead.tree.dens <- ifelse(dead.tree.dens >= 6, 1, 0) #assign MOG status (1) if condition meets forest type definition
                        
                        local.MOG.status <- ifelse(sum(stand.age,dead.tree.dens,live.large.tree.dens,stand.basal.area,na.rm=TRUE) == 4, 1, 0) #if all conditions are met, assign 1 for old growth. Otherwise assign 0 for non-OG
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end instructions for (southern) eastern riverfront forest
                      
                      if(southern.type[s] %in% 29){ #begin instructions for (southern) southern wet pine forest, woodland, and savanna
                        
                        #determine if information meets MOG criteria
                        stand.age <- ifelse(stand.age >= 80, 1, 0)
                        
                        large.trees <- live.trees[which(live.trees$DIA >= 9), ] #isolate large trees per the forest type definition of large
                        if(nrow(large.trees) > 0){live.large.tree.dens <- nrow(large.trees)/condition.area} #calculate density of large trees
                        if(nrow(large.trees) == 0){live.large.tree.dens <- 0} #handle large tree density if no large trees exist
                        live.large.tree.dens <- ifelse(live.large.tree.dens >= 6, 1, 0) #assign MOG status (1) if condition meets forest type definition
                        
                        stand.basal.area <- ifelse(stand.basal.area >= 10, 1, 0) #assign MOG status (1) if condition meets forest type definition
                        
                        dead.tree.dens <- ifelse(dead.tree.dens >= 0, 1, 0) #assign MOG status (1) if condition meets forest type definition
                        
                        local.MOG.status <- ifelse(sum(stand.age+dead.tree.dens,live.large.tree.dens,stand.basal.area,na.rm=TRUE) == 4, 1, 0) #if all conditions are met, assign 1 for old growth. Otherwise assign 0 for non-OG
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end instructions for (southern) southern wet pine forest, woodland, and savanna
                      
                      if(southern.type[s] %in% 31){ #begin instructions for (southern) montane and allied spruce and spruce-fir forest
                        
                        #determine if information meets MOG criteria
                        stand.age <- ifelse(stand.age >= 120, 1, 0)
                        
                        large.trees <- live.trees[which(live.trees$DIA >= 20), ] #isolate large trees per the forest type definition of large
                        if(nrow(large.trees) > 0){live.large.tree.dens <- nrow(large.trees)/condition.area} #calculate density of large trees
                        if(nrow(large.trees) == 0){live.large.tree.dens <- 0} #handle large tree density if no large trees exist
                        live.large.tree.dens <- ifelse(live.large.tree.dens >= 6, 1, 0) #assign MOG status (1) if condition meets forest type definition
                        
                        stand.basal.area <- ifelse(stand.basal.area >= 40, 1, 0) #assign MOG status (1) if condition meets forest type definition
                        
                        dead.tree.dens <- ifelse(dead.tree.dens >= 14, 1, 0) #assign MOG status (1) if condition meets forest type definition
                        
                        local.MOG.status <- ifelse(sum(stand.age,dead.tree.dens,live.large.tree.dens,stand.basal.area,na.rm=TRUE) == 4, 1, 0) #if all conditions are met, assign 1 for old growth. Otherwise assign 0 for non-OG
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end instructions for (southern) montane and allied spruce and spruce-fir forest
                      
                      #if old-growth is not designated, assess if the plot is at least mature
                      if(local.MOG.status == 0){ #begin mature assessment based on Table 19 of MOG document
                        
                        if(forest.type %in% c(105,162,402,407,401,406,409,405)){#begin maturity assessment for (southern) conifer southern hardwood
                          #assess of indices meet maturity standards for specific forest type
                          t.QMDdom <- ifelse(QMDdom >= 8.3, 0.42, 0)
                          t.tpadom <- ifelse(tpadom <= 111.6, 0.3, 0)
                          t.HTquart <- ifelse(HTquart >= 39.2, 0.28, 0)
                          
                          local.MOG.status <- sum(t.QMDdom,t.tpadom,t.HTquart,na.rm=TRUE) 
                          MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                          
                        }#end maturity assessment for (southern) conifer southern hardwood
                        
                        if(forest.type %in% c(141,403)){#begin maturity assessment for (southern) longleaf pine
                          #assess of indices meet maturity standards for specific forest type
                          t.QMDdom <- ifelse(QMDdom >= 10.2, 0.31, 0)
                          t.ddiscore <- ifelse(ddiscore >= 19, 0.23, 0)
                          t.tpadom <- ifelse(tpadom <= 54.7, 0.23, 0)
                          t.HTsd <- ifelse(HTsd >= 24, 0.12, 0)
                          t.badom <- ifelse(badom >= 44.7, 0.12, 0)
                          
                          local.MOG.status <- sum(t.QMDdom,t.ddiscore,t.tpadom,t.HTsd,t.badom,na.rm=TRUE)
                          MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                          
                        }#end maturity assessment for (southern) longleaf pine
                        
                        if(forest.type %in% c(502,510,515,514,505,504,503,501)){#begin maturity assessment for (southern) oaks
                          #assess of indices meet maturity standards for specific forest type
                          t.QMDdom <- ifelse(QMDdom >= 9.5, 0.3, 0)
                          t.ddiscore <- ifelse(ddiscore >= 22.8, 0.28, 0)
                          t.HTquart <- ifelse(HTquart >= 44.1, 0.22, 0)
                          t.badom <- ifelse(badom >= 55, 0.2, 0)
                          
                          local.MOG.status <- sum(t.QMDdom,t.ddiscore,t.HTquart,t.badom,na.rm=TRUE)
                          MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                          
                        }#end maturity assessment for (southern) oaks
                        
                        if(forest.type %in% c(103,104,166,142,123,165,161,164,163,167,162)){#begin maturity assessment for (southern) pines-conifers
                          #assess of indices meet maturity standards for specific forest type
                          t.QMDdom <- ifelse(QMDdom >= 11.4, 0.38, 0)
                          t.tpadom <- ifelse(tpadom <= 60.4, 0.26, 0)
                          t.HTquart <- ifelse(HTquart >= 65.8, 0.19, 0)
                          t.HTsd <- ifelse(HTsd >= 38.6, 0.17, 0)
                          
                          local.MOG.status <- sum(t.QMDdom,t.tpadom,t.HTquart,t.HTsd,na.rm=TRUE)
                          MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                          
                        }#end maturity assessment for (southern) pines-conifers
                        
                        if(forest.type %in% c(609,520,507,516,708,608,607,962,801,703,519,602,511,802,605,706,517,809,508,506,512,905,601,704,805,702,705)){#begin maturity assessment for (southern) southern hardwoods
                          #assess of indices meet maturity standards for specific forest type
                          t.ddiscore <- ifelse(ddiscore >= 30.1, 0.31, 0)
                          t.HTquart <- ifelse(HTquart >= 43.8, 0.26, 0)
                          t.badom <- ifelse(badom >= 59.1, 0.22, 0)
                          t.HTsd <- ifelse(HTsd >= 48, 0.21, 0)
                          
                          local.MOG.status <- sum(t.ddiscore,t.HTquart,t.badom,t.HTsd,na.rm=TRUE)
                          MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                          
                        }#end maturity assessment for (southern) southern hardwoods
                        
                      } #end mature assessment
                      
                    }#end looping through southern forest types (s)
                    
                  } #end instructions if condition fits within a defined southern forest type
                  
                }#end instructions for southern region
                
                if(region %in% "rocky"){#begin instructions for rock mountain region (based on Table 8 of MOG document)
                  
                  n.trees.broken.p.acre <- nrow(local.trees[local.trees$CULLMSTOP > 0 | local.trees$CULL > 0, ])/condition.area #calculate number of trees with cull or broken top
                  n.dead.trees <- nrow(local.trees[local.trees$DIA >= 10 & local.trees$STATUSCD == 2, ])/condition.area #calculate the number of dead trees with dbh of at least 10 inches
                  
                  #pull tree age from various aspects of the dataset
                  local.trees$tree.age <- ifelse(!is.na(local.trees$BHAGE),local.trees$BHAGE, 
                                                 ifelse(!is.na(local.trees$TOTAGE),local.trees$TOTAGE,raw.stand.age))
                  
                  if(forest.type %in% c(220,221,222,224,225,226,200,201,202,203,120:129,260:271)){#begin old growth assessment for (rocky mountain) ponderosa pine & mixed conifer &spruce/fir
                    
                    local.large.trees <- local.trees[local.trees$DIA >= 16 & local.trees$STATUSCD == 1, ] #isolate to large living tress only, using forest-type cutoff for dbh
                    local.large.old.trees <- local.large.trees[local.large.trees$tree.age >= 200, ]
                    
                    n.large.old.trees.p.acre <- nrow(local.large.old.trees)/condition.area #convert from count to n per acre
                    n.large.old.trees.p.acre <- ifelse(n.large.old.trees.p.acre >= 10, 1, 0) #assess if stand meets forest-type specific old growth threshold
                    
                    n.trees.broken.p.acre <- ifelse(n.trees.broken.p.acre >= 1, 1, 0) #assess if stand meets forest-type specific old growth threshold
                    n.dead.trees <- ifelse(n.dead.trees >= 2, 1, 0) #assess if stand meets forest-type specific old growth threshold
                    
                    local.MOG.status <- ifelse(sum(n.large.old.trees.p.acre,n.trees.broken.p.acre,n.dead.trees,na.rm=TRUE) == 3, 1, 0) #if the plot meets all old-growth thresholds, assign 1 for old growth. If not, assign 0 for not old growth.
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (rocky mountain) ponderosa pine & mixed conifer &spruce/fir
                  
                  if(forest.type %in% c(900:905)){#begin old growth assessment for (rocky mountain) aspen
                    
                    local.large.trees <- local.trees[local.trees$DIA >= 14 & local.trees$STATUSCD == 1, ] #isolate to large living tress only, using forest-type cutoff for dbh
                    local.large.old.trees <- local.large.trees[local.large.trees$tree.age >= 200, ]
                    
                    n.large.old.trees.p.acre <- nrow(local.large.old.trees)/condition.area #convert from count to n per acre
                    n.large.old.trees.p.acre <- ifelse(n.large.old.trees.p.acre >= 10, 1, 0) #assess if stand meets forest-type specific old growth threshold
                    
                    n.trees.broken.p.acre <- ifelse(n.trees.broken.p.acre >= 1, 1, 0) #assess if stand meets forest-type specific old growth threshold
                    n.dead.trees <- ifelse(n.dead.trees >= 0, 1, 0) #assess if stand meets forest-type specific old growth threshold
                    
                    local.MOG.status <- ifelse(sum(n.large.old.trees.p.acre,n.trees.broken.p.acre,n.dead.trees,na.rm=TRUE) == 3, 1, 0) #if the plot meets all old-growth thresholds, assign 1 for old growth. If not, assign 0 for not old growth.
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (rocky mountain) aspen
                  
                  if(forest.type %in% c(280,281)){#begin old growth assessment for (rocky mountain) lodgepole pine
                    
                    local.large.trees <- local.trees[local.trees$DIA >= 10 & local.trees$STATUSCD == 1, ] #isolate to large living tress only, using forest-type cutoff for dbh
                    local.large.old.trees <- local.large.trees[local.large.trees$tree.age >= 150, ]
                    
                    n.large.old.trees.p.acre <- nrow(local.large.old.trees)/condition.area #convert from count to n per acre
                    n.large.old.trees.p.acre <- ifelse(n.large.old.trees.p.acre >= 10, 1, 0) #assess if stand meets forest-type specific old growth threshold
                    
                    n.trees.broken.p.acre <- ifelse(n.trees.broken.p.acre >= 1, 1, 0) #assess if stand meets forest-type specific old growth threshold
                    n.dead.trees <- ifelse(n.dead.trees >= 2, 1, 0) #assess if stand meets forest-type specific old growth threshold
                    
                    local.MOG.status <- ifelse(sum(n.large.old.trees.p.acre,n.trees.broken.p.acre,n.dead.trees,na.rm=TRUE) == 3, 1, 0) #if the plot meets all old-growth thresholds, assign 1 for old growth. If not, assign 0 for not old growth.
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (rocky mountain) lodgepole pine
                  
                  if(forest.type %in% c(180,182,184,185)){#begin old growth assessment for (rocky mountain) pinyon-juniper
                    
                    local.large.trees <- local.trees[local.trees$DIA >= 12 & local.trees$STATUSCD == 1, ] #isolate to large living tress only, using forest-type cutoff for dbh
                    local.large.old.trees <- local.large.trees[local.large.trees$tree.age >= 200, ]
                    
                    n.large.old.trees.p.acre <- nrow(local.large.old.trees)/condition.area #convert from count to n per acre
                    n.large.old.trees.p.acre <- ifelse(n.large.old.trees.p.acre >= 30, 1, 0) #assess if stand meets forest-type specific old growth threshold
                    
                    n.trees.broken.p.acre <- ifelse(n.trees.broken.p.acre >= 1, 1, 0) #assess if stand meets forest-type specific old growth threshold
                    n.dead.trees <- ifelse(n.dead.trees >= 1, 1, 0) #assess if stand meets forest-type specific old growth threshold
                    
                    local.MOG.status <- ifelse(sum(n.large.old.trees.p.acre,n.trees.broken.p.acre,n.dead.trees,na.rm=TRUE) == 3, 1, 0) #if the plot meets all old-growth thresholds, assign 1 for old growth. If not, assign 0 for not old growth.
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (rocky mountain) pinyon-juniper
                  
                  if(forest.type %in% c(360:369)){#begin old growth assessment for (rocky mountain) white pine
                    
                    local.large.trees <- local.trees[local.trees$DIA >= 12 & local.trees$STATUSCD == 1, ] #isolate to large living tress only, using forest-type cutoff for dbh
                    local.large.old.trees <- local.large.trees[local.large.trees$tree.age >= 200, ]
                    
                    n.large.old.trees.p.acre <- nrow(local.large.old.trees)/condition.area #convert from count to n per acre
                    n.large.old.trees.p.acre <- ifelse(n.large.old.trees.p.acre >= 10, 1, 0) #assess if stand meets forest-type specific old growth threshold
                    
                    n.trees.broken.p.acre <- ifelse(n.trees.broken.p.acre >= 0, 1, 0) #assess if stand meets forest-type specific old growth threshold
                    n.dead.trees <- ifelse(n.dead.trees >= 0, 1, 0) #assess if stand meets forest-type specific old growth threshold
                    
                    local.MOG.status <- ifelse(sum(n.large.old.trees.p.acre,n.trees.broken.p.acre,n.dead.trees,na.rm=TRUE) == 3, 1, 0) #if the plot meets all old-growth thresholds, assign 1 for old growth. If not, assign 0 for not old growth.
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (rocky mountain) white pine
                  
                  if(forest.type %in% c(970:976)){#begin old growth assessment for (rocky mountain) gambel oak
                    
                    local.large.trees <- local.trees[local.trees$DIA >= 4 & local.trees$STATUSCD == 1, ] #isolate to large living tress only, using forest-type cutoff for dbh
                    local.large.old.trees <- local.large.trees[local.large.trees$tree.age >= 80, ]
                    
                    n.large.old.trees.p.acre <- nrow(local.large.old.trees)/condition.area #convert from count to n per acre
                    n.large.old.trees.p.acre <- ifelse(n.large.old.trees.p.acre >= 30, 1, 0) #assess if stand meets forest-type specific old growth threshold
                    
                    n.trees.broken.p.acre <- ifelse(n.trees.broken.p.acre >= 0, 1, 0) #assess if stand meets forest-type specific old growth threshold
                    n.dead.trees <- ifelse(n.dead.trees >= 0, 1, 0) #assess if stand meets forest-type specific old growth threshold
                    
                    local.MOG.status <- ifelse(sum(n.large.old.trees.p.acre,n.trees.broken.p.acre,n.dead.trees,na.rm=TRUE) == 3, 1, 0) #if the plot meets all old-growth thresholds, assign 1 for old growth. If not, assign 0 for not old growth.
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (rocky mountain) gambel oak
                  
                  if(forest.type %in% c(700:722)){#begin old growth assessment for (rocky mountain) cottonwood
                    
                    local.large.trees <- local.trees[local.trees$DIA >= 14 & local.trees$STATUSCD == 1, ] #isolate to large living tress only, using forest-type cutoff for dbh
                    local.large.old.trees <- local.large.trees[local.large.trees$tree.age >= 100, ]
                    
                    n.large.old.trees.p.acre <- nrow(local.large.old.trees)/condition.area #convert from count to n per acre
                    n.large.old.trees.p.acre <- ifelse(n.large.old.trees.p.acre >= 20, 1, 0) #assess if stand meets forest-type specific old growth threshold
                    
                    n.trees.broken.p.acre <- ifelse(n.trees.broken.p.acre >= 0, 1, 0) #assess if stand meets forest-type specific old growth threshold
                    n.dead.trees <- ifelse(n.dead.trees >= 0, 1, 0) #assess if stand meets forest-type specific old growth threshold
                    
                    local.MOG.status <- ifelse(sum(n.large.old.trees.p.acre,n.trees.broken.p.acre,n.dead.trees,na.rm=TRUE) == 3, 1, 0) #if the plot meets all old-growth thresholds, assign 1 for old growth. If not, assign 0 for not old growth.
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (rocky mountain) cottonwood
                  
                  #if not old growth, check to see if the site is at least mature (based on Table 19 of MOG document)
                  if(local.MOG.status == 0){ #begin (rocky mountain) maturity assessment
                    
                    if(forest.type %in% c(900:905,703,500:520,961,962)){#start maturity assessment for (rocky mountain) aspen/cottonwood/oaks
                      
                      t.HTquart <- ifelse(HTquart >= 32.9, 0.31, 0)
                      t.ddiscore <- ifelse(ddiscore >= 18.6, 0.27, 0)
                      t.badom <- ifelse(badom >= 55.1, 0.26, 0)
                      t.HTsd <- ifelse(HTsd >= 25.3, 0.15, 0)
                      
                      local.MOG.status <- sum(t.HTquart,t.ddiscore,t.badom,t.HTsd, na.rm = TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end maturity assessment for (rocky mountain) aspen/cottonwood/oaks
                    
                    
                    if(forest.type %in% c(201)){#start maturity assessment for (rocky mountain) douglas fir
                      
                      t.ddiscore <- ifelse(ddiscore >= 29.2, 0.3, 0)
                      t.badom <- ifelse(badom >= 65.8, 0.21, 0)
                      t.HTquart <- ifelse(HTquart >= 40.6, 0.18, 0)
                      t.QMDdom <- ifelse(QMDdom >= 9.3, 0.17, 0)
                      t.snagbatot <- ifelse(snagbatot >= 21.3, 0.15, 0)
                      
                      local.MOG.status <- sum(t.ddiscore,t.badom,t.HTquart,t.QMDdom,t.snagbatot,na.rm=TRUE) 
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end maturity assessment for (rocky mountain) douglas fir
                    
                    if(forest.type %in% c(970:976)){#start maturity assessment for (rocky mountain) gambel oak
                      
                      badom <- ifelse(badom >= 25.3, 0.3, 0)
                      ddiscore <- ifelse(ddiscore >= 8, 0.25, 0)
                      HTquart <- ifelse(HTquart >= 10.4, 0.24, 0)
                      QMDdom <- ifelse(QMDdom >= 2.9, 0.21, 0)
                      
                      local.MOG.status <- ifelse(sum(badom,ddiscore,HTquart,QMDdom, na.rm=TRUE) >= 0.5, 0.5, 0) #if score is greater than 0.5, assign 0.5 for mature. If not, assign 0 for not mature or old growth
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end maturity assessment for (rocky mountain) gambel oak
                    
                    if(forest.type %in% c(280,281)){#start maturity assessment for (rocky mountain) lodgepole pine
                      
                      t.QMDdom <- ifelse(QMDdom >= 3.7, 0.56, 0)
                      t.badom <- ifelse(badom >= 33.8, 0.38, 0)
                      t.HTsd <- ifelse(HTsd >= 17.5, 0.16, 0)
                      
                      local.MOG.status <- sum(t.QMDdom,t.badom,t.HTsd,na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end maturity assessment for (rocky mountain) lodgepole pine
                    
                    if(forest.type %in% c(360:369,170,171,172)){#start maturity assessment for (rocky mountain) other western softwoods
                      
                      t.ddiscore <- ifelse(ddiscore >= 24, 0.32, 0)
                      t.QMDdom <- ifelse(QMDdom >= 6.5, 0.29, 0)
                      t.HTquart <- ifelse(HTquart >= 28.2, 0.24, 0)
                      t.HTsd <- ifelse(HTsd >= 21.6, 0.15, 0)
                      
                      local.MOG.status <- sum(t.ddiscore,t.QMDdom,t.HTquart,t.HTsd, na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end maturity assessment for (rocky mountain) other western softwoods
                    
                    if(forest.type %in% c(180,182,184,185)){#start maturity assessment for (rocky mountain) pinyon juniper
                      
                      t.ddiscore <- ifelse(ddiscore >= 33.5, 0.55, 0)
                      t.QMDdom <- ifelse(QMDdom >= 8.6, 0.45, 0)
                      
                      local.MOG.status <- sum(t.ddiscore,t.QMDdom, na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end maturity assessment for (rocky mountain) other pinyon juniper
                    
                    if(forest.type %in% c(220,221,222,224,225,226)){#start maturity assessment for (rocky mountain) ponderosa pine
                      
                      t.QMDdom <- ifelse(QMDdom >= 11.8, 0.33, 0)
                      t.ddiscore <- ifelse(ddiscore >= 31.6, 0.28, 0)
                      t.HTsd <- ifelse(HTsd >= 39, 0.21, 0)
                      t.badom <- ifelse(badom >= 67.3, 0.18, 0)
                      
                      local.MOG.status <- sum(t.QMDdom,t.ddiscore,t.HTsd,t.badom, na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end maturity assessment for (rocky mountain) ponderosa pine
                    
                    if(forest.type %in% c(260:271,120:129)){#start maturity assessment for (rocky mountain) spruce/fir
                      
                      t.ddiscore <- ifelse(ddiscore >= 28.8, 0.31, 0)
                      t.badom <- ifelse(badom >= 87.2, 0.27, 0)
                      t.HTquart <- ifelse(HTquart >= 43.5, 0.24, 0)
                      t.HTsd <- ifelse(HTsd >= 44.6, 0.18, 0)
                      
                      local.MOG.status <- sum(t.ddiscore,t.badom,t.HTquart,t.HTsd, na.rm = TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end maturity assessment for (rocky mountain) spruce/fir
                  } #end (rocky mountain) maturity assessment
                  
                  #conditions[c,"MOG"] <- max(MOG.vector,na.rm=TRUE) #if MOG conditions are met for any of the forest types assigned to the plot, assign 1 for OG. Otherwise assign 0.5 for mature or zero for neither
                  #conditions[c,"age"] <- raw.stand.age #save stand age
                  
                }#end instructions for rocky mountain region
                
                if(region %in% "intermountain"){#begin instructions for intermountain region
                  
                  #assign intermountain old growth type (based on Table 11 of MOG document)
                  intermountain.type <- NA #assign type as NA to trigger additional data download downstream if none of the following conditions hold
                  if(!(local.plot$STATECD %in% 16) & local.condition$FORTYPCD %in% c(265,261)){intermountain.type <- "engelmann spruce - subalpine fir - warm - UT"}
                  if(!(local.plot$STATECD %in% 16) & local.condition$FORTYPCD %in% 266 & local.condition$PHYSCLCD > 20){intermountain.type <- "engelmann spruce - subalpine fir - warm - UT"}
                  if(!(local.plot$STATECD %in% 16) & local.condition$FORTYPCD %in% 266 & !any(local.trees$SPCD %in% c(113,101,72))){intermountain.type <- "engelmann spruce - subalpine fir - warm - UT"}
                  if(local.plot$STATECD %in% 16 & local.condition$FORTYPCD %in% c(265,261)){intermountain.type <- "engelmann spruce - subalpine fir - warm - ID"}
                  if(local.plot$STATECD %in% 16 & local.condition$FORTYPCD %in% 266 & local.condition$PHYSCLCD > 20){intermountain.type <- "engelmann spruce - subalpine fir - warm - ID"}
                  if(local.plot$STATECD %in% 16 & local.condition$FORTYPCD %in% 266 & !any(local.trees$SPCD %in% c(113,101,72))){intermountain.type <- "engelmann spruce - subalpine fir - warm - ID"}
                  if(local.condition$FORTYPCD %in% 266 & local.condition$PHYSCLCD < 20){intermountain.type <- "engelmann spruce - subalpine fir - cold"}
                  if(local.condition$FORTYPCD %in% 266 & any(local.trees$SPCD %in% c(113,101,72))){intermountain.type <- "engelmann spruce - subalpine fir - cold"}
                  if(local.condition$FORTYPCD %in% 268 & local.condition$SITECLCD < 7){intermountain.type <- "engelmann spruce - subalpine fir - cold"}
                  if(local.condition$FORTYPCD %in% 268 & local.condition$SITECLCD == 7){intermountain.type <- "engelmann spruce - subalpine fir - alpine"}
                  if(local.condition$FORTYPCD %in% 367){intermountain.type <- "whitebark pine"}
                  if(local.condition$FORTYPCD %in% 365){intermountain.type <- "bristlecone pine"}
                  if(local.condition$FORTYPCD %in% 201 & local.condition$SITECLCD < 6){intermountain.type <- "douglas - fir - high"}
                  if(local.condition$FORTYPCD %in% 201 & local.condition$SITECLCD >= 6){intermountain.type <- "douglas - fir - low"}
                  if(local.condition$FORTYPCD %in% 267){intermountain.type <- "grand fir"}
                  if(local.condition$FORTYPCD %in% 269){intermountain.type <- "blue spruce"}
                  if(local.condition$FORTYPCD %in% c(371,262) & local.condition$PHYSCLCD < 20){intermountain.type <- "conifer mixed forest - low"}
                  if(local.condition$FORTYPCD %in% c(371,262) & local.condition$PHYSCLCD >= 20){intermountain.type <- "conifer mixed forest - productive"}
                  if(local.condition$FORTYPCD %in% 901 & local.condition$PHYSCLCD < 20){intermountain.type <- "aspen - dry"}
                  if(local.condition$FORTYPCD %in% 901 & local.condition$PHYSCLCD > 20){intermountain.type <- "aspen - mesic"}
                  if(local.condition$FORTYPCD %in% 281){intermountain.type <- "lodgepole pine"}
                  if(local.condition$FORTYPCD %in% 366 & local.condition$SITECLCD > 6){intermountain.type <- "limber pine - lower"}
                  if(local.condition$FORTYPCD %in% 366 & local.condition$SITECLCD <= 6){intermountain.type <- "limber pine - montane"}
                  if(local.condition$FORTYPCD %in% c(220,221,222,225) & local.condition$ADFORCD %in% c(402,412,413,414) & local.condition$SITECLCD > 5){intermountain.type <- "ponderosa pine - n - seral"}
                  if(local.condition$FORTYPCD %in% c(220,221,222,225) & local.condition$ADFORCD %in% c(402,412,413,414) & local.condition$SITECLCD <= 5){intermountain.type <- "ponderosa pine - n - climax"}
                  if(local.condition$FORTYPCD %in% c(220,221,222,225) & !(local.condition$ADFORCD %in% c(402,412,413,414)) & local.condition$SITECLCD > 5){intermountain.type <- "ponderosa pine - rm - seral"}
                  if(local.condition$FORTYPCD %in% c(220,221,222,225) & !(local.condition$ADFORCD %in% c(402,412,413,414)) & local.condition$SITECLCD <= 5){intermountain.type <- "ponderosa pine - rm - climax"}
                  if(local.condition$FORTYPCD %in% c(182,184,185) & local.condition$ADFORCD %in% c(402,403,412,413,414,415,417,420) & local.condition$PHYSCLCD < 20){intermountain.type <- "pinyon - juniper - nw - low"}
                  if(local.condition$FORTYPCD %in% c(182,184,185) & local.condition$ADFORCD %in% c(402,403,412,413,414,415,417,420) & local.condition$PHYSCLCD > 20){intermountain.type <- "pinyon - juniper - nw - high"}
                  if(local.condition$FORTYPCD %in% c(182,184,185) & local.condition$ADFORCD %in% c(401,407,408,410) & local.condition$PHYSCLCD < 20){intermountain.type <- "pinyon - juniper - se - low"}
                  if(local.condition$FORTYPCD %in% c(182,184,185) & local.condition$ADFORCD %in% c(401,407,408,410) & local.condition$PHYSCLCD > 20){intermountain.type <- "pinyon - juniper - se - high"}
                  
                  if(is.na(intermountain.type)){#begin instructions if an intermountain forest type has not yet been assigned
                    
                    ECOSUBCD <- FIA.geomplot[which(FIA.geomplot$CN %in% local.plot$CN),"ECOSUBCD"] #assign ecological subregion
                    
                    if(local.condition$FORTYPCD %in% c(182,184,185) & local.condition$ADFORCD %in% c(418,419) & ECOSUBCD %in% c("M331Dn","M331Do","M331Dv","M331Di") & local.condition$PHYSCLCD < 20){intermountain.type <- "pinyon - juniper - nw - low"}
                    if(local.condition$FORTYPCD %in% c(182,184,185) & local.condition$ADFORCD %in% c(418,419) & ECOSUBCD %in% c("M331Dn","M331Do","M331Dv","M331Di") & local.condition$PHYSCLCD > 20){intermountain.type <- "pinyon - juniper - nw - high"}
                    if(local.condition$FORTYPCD %in% c(182,184,185) & local.condition$ADFORCD %in% c(418,419) & !(ECOSUBCD %in% c("M331Dn","M331Do","M331Dv","M331Di")) & local.condition$PHYSCLCD < 20){intermountain.type <- "pinyon - juniper - se - low"}
                    if(local.condition$FORTYPCD %in% c(182,184,185) & local.condition$ADFORCD %in% c(418,419) & !(ECOSUBCD %in% c("M331Dn","M331Do","M331Dv","M331Di")) & local.condition$PHYSCLCD > 20){intermountain.type <- "pinyon - juniper - se - high"}
                    
                  }#end instructions if an intermountain forest type has not yet been assigned
                  
                  #pull tree age from various aspects of the dataset
                  local.trees$tree.age <- ifelse(!is.na(local.trees$BHAGE),local.trees$BHAGE, 
                                                 ifelse(!is.na(local.trees$TOTAGE),local.trees$TOTAGE,raw.stand.age)) #if no data to the contrary, assume trees are the age of the stand itself
                  
                  
                  #asses old growth conditions based on Table 11 of MOG document
                  if(intermountain.type %in% "engelmann spruce - subalpine fir - warm - UT"){#begin old growth assessment for (intermountain) "engelmann spruce - subalpine fir - warm - UT"
                    
                    local.large.old.trees <- local.trees[local.trees$tree.age >= 220 & local.trees$DIA >= 20, ]
                    local.large.old.trees.p.acre <- nrow(local.large.old.trees)/condition.area
                    
                    local.MOG.status <- ifelse(local.large.old.trees.p.acre >= 25, 1, 0) #assign old growth status (1) if threshold is met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (intermountain) "engelmann spruce - subalpine fir - warm - UT"
                  
                  if(intermountain.type %in% "engelmann spruce - subalpine fir - warm - ID"){#begin old growth assessment for (intermountain) "engelmann spruce - subalpine fir - warm - ID"
                    
                    local.large.old.trees <- local.trees[local.trees$tree.age >= 220 & local.trees$DIA >= 24, ]
                    local.large.old.trees.p.acre <- nrow(local.large.old.trees)/condition.area
                    
                    local.MOG.status <- ifelse(local.large.old.trees.p.acre >= 25, 1, 0) #assign old growth status (1) if threshold is met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (intermountain) "engelmann spruce - subalpine fir - warm - ID"
                  
                  if(intermountain.type %in% "engelmann spruce - subalpine fir - cold"){#begin old growth assessment for (intermountain) "engelmann spruce - subalpine fir - cold"
                    
                    local.large.old.trees <- local.trees[local.trees$tree.age >= 150 & local.trees$DIA >= 15, ]
                    local.large.old.trees.p.acre <- nrow(local.large.old.trees)/condition.area
                    
                    local.MOG.status <- ifelse(local.large.old.trees.p.acre >= 15, 1, 0) #assign old growth status (1) if threshold is met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (intermountain) "engelmann spruce - subalpine fir - cold"
                  
                  if(intermountain.type %in% "engelmann spruce - subalpine fir - alpine"){#begin old growth assessment for (intermountain) "engelmann spruce - subalpine fir - alpine"
                    
                    local.large.old.trees <- local.trees[local.trees$tree.age >= 150 & local.trees$DIA >= 12, ]
                    local.large.old.trees.p.acre <- nrow(local.large.old.trees)/condition.area
                    
                    local.MOG.status <- ifelse(local.large.old.trees.p.acre >= 10, 1, 0) #assign old growth status (1) if threshold is met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (intermountain) "engelmann spruce - subalpine fir - alpine"
                  
                  if(intermountain.type %in% "whitebark pine"){#begin old growth assessment for (intermountain) "whitebark pine"
                    
                    local.large.old.trees <- local.trees[local.trees$tree.age >= 250 & local.trees$DIA >= 18, ]
                    local.large.old.trees.p.acre <- nrow(local.large.old.trees)/condition.area
                    
                    local.MOG.status <- ifelse(local.large.old.trees.p.acre >= 15, 1, 0) #assign old growth status (1) if threshold is met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (intermountain) "whitebark pine"
                  
                  if(intermountain.type %in% "bristlecone pine"){#begin old growth assessment for (intermountain) "bristlecone pine"
                    
                    local.large.old.trees <- local.trees[local.trees$tree.age >= 300 & local.trees$DIA >= 10, ]
                    local.large.old.trees.p.acre <- nrow(local.large.old.trees)/condition.area
                    
                    local.MOG.status <- ifelse(local.large.old.trees.p.acre >= 5, 1, 0) #assign old growth status (1) if threshold is met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (intermountain) "bristlecone pine"
                  
                  if(intermountain.type %in% "douglas - fir - high"){#begin old growth assessment for (intermountain) "douglas - fir - high"
                    
                    local.large.old.trees <- local.trees[local.trees$tree.age >= 200 & local.trees$DIA >= 24, ]
                    local.large.old.trees.p.acre <- nrow(local.large.old.trees)/condition.area
                    
                    local.MOG.status <- ifelse(local.large.old.trees.p.acre >= 15, 1, 0) #assign old growth status (1) if threshold is met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (intermountain) "douglas - fir - high"
                  
                  if(intermountain.type %in% "douglas - fir - low"){#begin old growth assessment for (intermountain) "douglas - fir - low"
                    
                    local.large.old.trees <- local.trees[local.trees$tree.age >= 200 & local.trees$DIA >= 18, ]
                    local.large.old.trees.p.acre <- nrow(local.large.old.trees)/condition.area
                    
                    local.MOG.status <- ifelse(local.large.old.trees.p.acre >= 10, 1, 0) #assign old growth status (1) if threshold is met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (intermountain) "douglas - fir - low"
                  
                  if(intermountain.type %in% "grand fir"){#begin old growth assessment for (intermountain) "grand fir"
                    
                    local.large.old.trees <- local.trees[local.trees$tree.age >= 200 & local.trees$DIA >= 24, ]
                    local.large.old.trees.p.acre <- nrow(local.large.old.trees)/condition.area
                    
                    local.MOG.status <- ifelse(local.large.old.trees.p.acre >= 15, 1, 0) #assign old growth status (1) if threshold is met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (intermountain) "grand fir"
                  
                  if(intermountain.type %in% "blue spruce"){#begin old growth assessment for (intermountain) "blue spruce"
                    
                    local.large.old.trees <- local.trees[local.trees$tree.age >= 250 & local.trees$DIA >= 16, ]
                    local.large.old.trees.p.acre <- nrow(local.large.old.trees)/condition.area
                    
                    local.MOG.status <- ifelse(local.large.old.trees.p.acre >= 10, 1, 0) #assign old growth status (1) if threshold is met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (intermountain) "blue spruce"
                  
                  if(intermountain.type %in% "conifer mixed forest - low"){#begin old growth assessment for (intermountain) "conifer mixed forest - low"
                    
                    local.large.old.trees <- local.trees[local.trees$tree.age >= 256 & local.trees$DIA >= 29, ]
                    local.large.old.trees.p.acre <- nrow(local.large.old.trees)/condition.area
                    
                    local.MOG.status <- ifelse(local.large.old.trees.p.acre >= 11, 1, 0) #assign old growth status (1) if threshold is met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (intermountain) "conifer mixed forest - low"
                  
                  if(intermountain.type %in% "conifer mixed forest - productive"){#begin old growth assessment for (intermountain) "conifer mixed forest - productive"
                    
                    local.large.old.trees <- local.trees[local.trees$tree.age >= 188 & local.trees$DIA >= 39, ]
                    local.large.old.trees.p.acre <- nrow(local.large.old.trees)/condition.area
                    
                    local.MOG.status <- ifelse(local.large.old.trees.p.acre >= 10, 1, 0) #assign old growth status (1) if threshold is met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (intermountain) "conifer mixed forest - productive"
                  
                  if(intermountain.type %in% "aspen - dry"){#begin old growth assessment for (intermountain) "aspen - dry"
                    
                    local.large.old.trees <- local.trees[local.trees$tree.age >= 100 & local.trees$DIA >= 12, ]
                    local.large.old.trees.p.acre <- nrow(local.large.old.trees)/condition.area
                    
                    local.MOG.status <- ifelse(local.large.old.trees.p.acre >= 10, 1, 0) #assign old growth status (1) if threshold is met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (intermountain) "aspen - dry"
                  
                  if(intermountain.type %in% "aspen - mesic"){#begin old growth assessment for (intermountain) "aspen - mesic"
                    
                    local.large.old.trees <- local.trees[local.trees$tree.age >= 100 & local.trees$DIA >= 12, ]
                    local.large.old.trees.p.acre <- nrow(local.large.old.trees)/condition.area
                    
                    local.MOG.status <- ifelse(local.large.old.trees.p.acre >= 20, 1, 0) #assign old growth status (1) if threshold is met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (intermountain) "aspen - mesic"
                  
                  if(intermountain.type %in% "lodgepole pine"){#begin old growth assessment for (intermountain) "lodgepole pine"
                    
                    local.large.old.trees <- local.trees[local.trees$tree.age >= 140 & local.trees$DIA >= 11, ]
                    local.large.old.trees.p.acre <- nrow(local.large.old.trees)/condition.area
                    
                    local.MOG.status <- ifelse(local.large.old.trees.p.acre >= 25, 1, 0) #assign old growth status (1) if threshold is met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (intermountain) "lodgepole pine"
                  
                  if(intermountain.type %in% "limber pine - lower"){#begin old growth assessment for (intermountain) "limber pine - lower"
                    
                    local.large.old.trees <- local.trees[local.trees$tree.age >= 250 & local.trees$DIA >= 16, ]
                    local.large.old.trees.p.acre <- nrow(local.large.old.trees)/condition.area
                    
                    local.MOG.status <- ifelse(local.large.old.trees.p.acre >= 10, 1, 0) #assign old growth status (1) if threshold is met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (intermountain) "limber pine - lower"
                  
                  if(intermountain.type %in% "limber pine - montane"){#begin old growth assessment for (intermountain) "limber pine - montane"
                    
                    local.large.old.trees <- local.trees[local.trees$tree.age >= 500 & local.trees$DIA >= 16, ]
                    local.large.old.trees.p.acre <- nrow(local.large.old.trees)/condition.area
                    
                    local.MOG.status <- ifelse(local.large.old.trees.p.acre >= 10, 1, 0) #assign old growth status (1) if threshold is met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (intermountain) "limber pine - montane"
                  
                  if(intermountain.type %in% "ponderosa pine - n - seral"){#begin old growth assessment for (intermountain) "ponderosa pine - n - seral"
                    
                    local.large.old.trees <- local.trees[local.trees$tree.age >= 200 & local.trees$DIA >= 24, ]
                    local.large.old.trees.p.acre <- nrow(local.large.old.trees)/condition.area
                    
                    local.MOG.status <- ifelse(local.large.old.trees.p.acre >= 10, 1, 0) #assign old growth status (1) if threshold is met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (intermountain) "ponderosa pine - n - seral"
                  
                  if(intermountain.type %in% "ponderosa pine - n - climax"){#begin old growth assessment for (intermountain) "ponderosa pine - n - climax"
                    
                    local.large.old.trees <- local.trees[local.trees$tree.age >= 200 & local.trees$DIA >= 24, ]
                    local.large.old.trees.p.acre <- nrow(local.large.old.trees)/condition.area
                    
                    local.MOG.status <- ifelse(local.large.old.trees.p.acre >= 5, 1, 0) #assign old growth status (1) if threshold is met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (intermountain) "ponderosa pine - n - climax"
                  
                  if(intermountain.type %in% "ponderosa pine - rm - seral"){#begin old growth assessment for (intermountain) "ponderosa pine - rm - seral"
                    
                    local.large.old.trees <- local.trees[local.trees$tree.age >= 200 & local.trees$DIA >= 20, ]
                    local.large.old.trees.p.acre <- nrow(local.large.old.trees)/condition.area
                    
                    local.MOG.status <- ifelse(local.large.old.trees.p.acre >= 14, 1, 0) #assign old growth status (1) if threshold is met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (intermountain) "ponderosa pine - rm - seral"
                  
                  if(intermountain.type %in% "ponderosa pine - rm - climax"){#begin old growth assessment for (intermountain) "ponderosa pine - rm - climax"
                    
                    local.large.old.trees <- local.trees[local.trees$tree.age >= 200 & local.trees$DIA >= 16, ]
                    local.large.old.trees.p.acre <- nrow(local.large.old.trees)/condition.area
                    
                    local.MOG.status <- ifelse(local.large.old.trees.p.acre >= 7, 1, 0) #assign old growth status (1) if threshold is met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (intermountain) "ponderosa pine - rm - climax"
                  
                  if(intermountain.type %in% "pinyon - juniper - nw - low"){#begin old growth assessment for (intermountain) "pinyon - juniper - nw - low"
                    
                    local.large.old.trees <- local.trees[local.trees$tree.age >= 200 & local.trees$DIA >= 12, ]
                    local.large.old.trees.p.acre <- nrow(local.large.old.trees)/condition.area
                    
                    local.MOG.status <- ifelse(local.large.old.trees.p.acre >= 12, 1, 0) #assign old growth status (1) if threshold is met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (intermountain) "pinyon - juniper - nw - low"
                  
                  if(intermountain.type %in% "pinyon - juniper - nw - high"){#begin old growth assessment for (intermountain) "pinyon - juniper - nw - high"
                    
                    local.large.old.trees <- local.trees[local.trees$tree.age >= 250 & local.trees$DIA >= 18, ]
                    local.large.old.trees.p.acre <- nrow(local.large.old.trees)/condition.area
                    
                    local.MOG.status <- ifelse(local.large.old.trees.p.acre >= 30, 1, 0) #assign old growth status (1) if threshold is met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (intermountain) "pinyon - juniper - nw - high"
                  
                  if(intermountain.type %in% "pinyon - juniper - se - low"){#begin old growth assessment for (intermountain) "pinyon - juniper - se - low"
                    
                    local.large.old.trees <- local.trees[local.trees$tree.age >= 150 & local.trees$DIA >= 9, ]
                    local.large.old.trees.p.acre <- nrow(local.large.old.trees)/condition.area
                    
                    local.MOG.status <- ifelse(local.large.old.trees.p.acre >= 12, 1, 0) #assign old growth status (1) if threshold is met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (intermountain) "pinyon - juniper - se - low"
                  
                  if(intermountain.type %in% "pinyon - juniper - se - high"){#begin old growth assessment for (intermountain) "pinyon - juniper - se - high"
                    
                    local.large.old.trees <- local.trees[local.trees$tree.age >= 200 & local.trees$DIA >= 12, ]
                    local.large.old.trees.p.acre <- nrow(local.large.old.trees)/condition.area
                    
                    local.MOG.status <- ifelse(local.large.old.trees.p.acre >= 30, 1, 0) #assign old growth status (1) if threshold is met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (intermountain) "pinyon - juniper - se - high"
                  
                  
                  #If the plot is not old-growth, assess if it is at least mature based on Table 19 of MOG document
                  if(local.MOG.status == 0){#begin (intermountain) assessment of mature conditions
                    
                    if(intermountain.type %in% "aspen - dry"){#begin mature assessment for "aspen - dry" group
                      
                      t.badom <- ifelse(badom >= 22.6, 0.33, 0)
                      t.ddiscore <- ifelse(ddiscore >= 12, 0.3, 0)
                      t.HTquart <- ifelse(HTquart >= 16.8, 0.26, 0)
                      t.HTsd <- ifelse(HTsd >= 14.7, 0.11, 0)
                      
                      local.MOG.status <- sum(t.badom,t.ddiscore,t.HTquart,t.HTsd,na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for (intermountain) "aspen - dry" group
                    
                    if(intermountain.type %in% "aspen - mesic"){#begin mature assessment for "aspen - mesic" group
                      
                      t.HTquart <- ifelse(HTquart >= 29.1, 0.39, 0)
                      t.ddiscore <- ifelse(ddiscore >= 15.3, 0.37, 0)
                      t.HTsd <- ifelse(HTsd >= 16.3, 0.24, 0)
                      
                      local.MOG.status <- sum(t.HTquart,t.ddiscore,t.HTsd,na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for (intermountain) "aspen - mesic" group
                    
                    if(intermountain.type %in% c("bristlecone pine","limber pine - lower","limber pine - montane","whitebark pine")){#begin mature assessment for "bristlecone/limber/whitebark pines" group
                      
                      t.QMDdom <- ifelse(QMDdom >= 9.8, 0.39, 0)
                      t.badom <- ifelse(badom >= 77.2, 0.34, 0)
                      t.HTquart <- ifelse(HTquart >= 26.4, 0.27, 0)
                      
                      local.MOG.status <- sum(t.QMDdom,t.badom,t.HTquart,na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for (intermountain) "bristlecone/limber/whitebark pines" group
                    
                    if(intermountain.type %in% c("douglas - fir - high","douglas - fir - low","grand fir","conifer mixed forest - productive") | forest.type %in% c(260:271,320,321)){#begin mature assessment for "douglas fir" group
                      
                      t.ddiscore <- ifelse(ddiscore >= 33, 0.43, 0)
                      t.QMDdom <- ifelse(QMDdom >= 11.1, 0.34, 0)
                      t.HTquart <- ifelse(HTquart >= 40.2, 0.23, 0)
                      
                      local.MOG.status <- sum(t.ddiscore,t.QMDdom,t.HTquart,na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for (intermountain) "douglas fir" group
                    
                    if(forest.type %in% c(700:722)){#begin mature assessment for "elm/ash/cottonwood" group
                      
                      t.badom <- ifelse(badom >= 47.5, 0.42, 0)
                      t.ddiscore <- ifelse(ddiscore >= 19.1, 0.39, 0)
                      t.HTsd <- ifelse(HTsd >= 15.5, 0.18, 0)
                      
                      local.MOG.status <- sum(t.badom,t.ddiscore,t.HTsd,na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for (intermountain) "elm/ash/cottonwood" group
                    
                    if(intermountain.type %in% c("engelmann spruce - subalpine fir - warm - UT","engelmann spruce - subalpine fir - warm - ID","engelmann spruce - subalpine fir - alpine","blue spruce","engelmann spruce - subalpine fir - cold","conifer mixed forest - low")){#begin mature assessment for "engelmann spruce" group
                      
                      t.ddiscore <- ifelse(ddiscore >= 29.8, 0.32, 0)
                      t.QMDdom <- ifelse(QMDdom >= 8.3, 0.27, 0)
                      t.HTquart <- ifelse(HTquart >= 35.4, 0.23, 0)
                      t.HTsd <- ifelse(HTquart >= 57.4, 0.18, 0)
                      
                      local.MOG.status <- sum(t.ddiscore,t.QMDdom,t.HTquart,t.HTsd,na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for (intermountain) "engelmann spruce" group
                    
                    if(intermountain.type %in% "lodgepole pine"){#begin mature assessment for "lodgepole pine" group
                      
                      t.ddiscore <- ifelse(ddiscore >= 14.7, 0.3, 0)
                      t.HTquart <- ifelse(HTquart >= 23.5, 0.26, 0)
                      t.badom <- ifelse(badom >= 41.8, 0.26, 0)
                      t.HTsd <- ifelse(HTquart >= 18.1, 0.18, 0)
                      
                      local.MOG.status <- sum(t.ddiscore,t.HTquart,t.badom,t.HTsd,na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for (intermountain) "lodgepole pine" group
                    
                    if(intermountain.type %in% c("pinyon - juniper - nw - low","pinyon - juniper - nw - high") | forest.type %in% c(970:976,960:962,360:369)){#begin mature assessment for "pinyon juniper NW and others" group
                      
                      t.ddiscore <- ifelse(ddiscore >= 24, 0.42, 0)
                      t.QMDdom <- ifelse(QMDdom >= 8, 0.39, 0)
                      t.tpadom <- ifelse(tpadom <= 90.3, 0.19, 0)
                      
                      local.MOG.status <- sum(t.ddiscore,t.QMDdom,t.tpadom,na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for (intermountain) "pinyon juniper NW and others" group
                    
                    if(intermountain.type %in% "pinyon - juniper - se - high"){#begin mature assessment for "pinyon juniper SE high" group
                      
                      t.QMDdom <- ifelse(QMDdom >= 9.2, 0.52, 0)
                      t.ddiscore <- ifelse(ddiscore >= 32.9, 0.48, 0)
                      
                      local.MOG.status <- sum(t.QMDdom,t.ddiscore,na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for (intermountain) "pinyon juniper SE high" group
                    
                    if(intermountain.type %in% "pinyon - juniper - se - low"){#begin mature assessment for "pinyon juniper SE low" group
                      
                      t.QMDdom <- ifelse(QMDdom >= 8.3, 0.44, 0)
                      t.ddiscore <- ifelse(ddiscore >= 24, 0.56, 0)
                      
                      local.MOG.status <- sum(t.QMDdom,t.ddiscore,na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for (intermountain) "pinyon juniper SE low" group
                    
                    if(intermountain.type %in% c("ponderosa pine - n - climax","ponderosa pine - n - seral","ponderosa pine - rm - climax","ponderosa pine - rm - seral")){#begin mature assessment for "ponderosa pine" group
                      
                      t.QMDdom <- ifelse(QMDdom >= 14.2, 0.38, 0)
                      t.ddiscore <- ifelse(ddiscore >= 30.7, 0.22, 0)
                      t.HTquart <- ifelse(HTquart >= 49, 0.21, 0)
                      t.HTsd <- ifelse(HTsd >= 50.2, 0.19, 0)
                      
                      local.MOG.status <- sum(t.QMDdom,t.ddiscore,t.HTquart,t.HTsd,na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for (intermountain) "ponderosa pine" group
                  }#end (intermountain) assessment of mature conditions
                  
                }#end instructions for intermountain region
                
                if(region %in% "northern"){#begin instructions for northern region
                  
                  #determine dominant tree
                  tree.spp <- as.vector(unique(local.trees$SPCD)) #id
                  dom.trees <- data.frame() #create empty data frame to fill
                  keep.veg <- c("AGSP","FEID","FESC","PURT","SYAL","BERE","PRIV","SHCA","PHMA","SMST","CARU","VAGL","ARUV","XETE","VACA",
                                "CAGE","SPBE","JUCO","ARCO","SYOR","CLUN","ASCA","MEFE","ARNU","SETR","TABR","LIBO","OPHO","ATFI","ADPE",
                                "SYAL","COOC","GYDR","EQAR","GATR","STAM","LICA","CACA","LEGL","LUHI","VASC","CLPS","PSME","RIMO","ABLA",
                                "FIED","PUTR","PRVI","ASCR","SEST","ALSI","ALFI","UBO","THOC","GAGE","JUSC","MARE","SYOC","PSSP","PIPO",
                                "PIAB","PIBR","PIEN","PIGL","PIMA","PIPU","PIRU","PISI","PIGL","PIOM")
                  
                  for(q in 1:length(tree.spp)){#begin looping through trees
                    
                    local.subset <- local.trees[local.trees$SPCD %in% tree.spp[q], ] #subset to species of interest
                    
                    #add data to data frame
                    dom.trees[q,"SPCD"] <- tree.spp[q]
                    dom.trees[q,"dominance"] <- mean(local.subset$CCLCD, na.rm=TRUE) #calculate mean crown class code
                    
                  } #end looping through trees
                  
                  #join dominant tree data frame to list of trees to obtain species code
                  dom.trees <- merge(x = dom.trees,
                                     y = USFStrees,
                                     by.x = "SPCD",
                                     by.y = "FIA.Code",
                                     all.x = TRUE, all.y = FALSE)
                  veg.code <- dom.trees[which(dom.trees$dominance %in% min(dom.trees$dominance, na.rm=TRUE)), "PLANTS.Code"] #isolate species code for dominant tree
                  veg.code <- substr(veg.code,0,4) #retain only the first four characters of code (gets rid of trailing numbers)
                  
                  if(length(veg.code)>1){veg.code <- veg.code[1]} #if there are two equally dominant trees, arbitrarily select the first 
                  
                  #determine dominant understory
                  local.veg <- FIA.veg[FIA.veg$PLT_CN %in% local.plot$CN, ] #subset vegetation data to the plot level
                  local.veg <- local.veg[local.veg$CONDID %in% local.condition$CONDID, ] #subset vegetation data to the condition level
                  local.veg <- local.veg[(substr(local.veg$VEG_SPCD,0,4) %in% keep.veg), ] #subset to only understory plants that are used to classify habitat groups
                  veg.spp <- as.vector(unique(local.veg$VEG_SPCD)) #id
                  veg.spp <- substr(veg.spp,0,4) #remove any numbers appended to vegetation
                  dom.veg <- data.frame() #create empty data frame to fill
                  
                  if(length(veg.spp)>0){ #begin instructions for if vegetation data exists
                    for(q in 1:length(veg.spp)){#begin looping through trees
                      
                      local.subset <- local.veg[local.veg$VEG_SPCD %in% veg.spp[q], ] #subset to species of interest
                      
                      #add data to data frame
                      dom.veg[q,"SPCD"] <- veg.spp[q]
                      dom.veg[q,"dominance"] <- mean(local.subset$COVER_PCT, na.rm=TRUE) #calculate mean crown class code
                      
                    } #end looping through veg
                    
                    dominant.veg <- dom.veg[which(dom.veg$dominance %in% max(dom.veg$dominance)), "SPCD"] #isolate species code for dominant understory
                    if(length(dominant.veg)==0){dominant.veg <- NA} #if nothing is dominant, save as NA
                    if(length(dominant.veg)>1){dominant.veg <- dominant.veg[1]} #if multiple under story species are co-dominant, arbitrarily select the first (classifications are expanded below to capture all combinations, so arbitrary selection does not impact classificatoin process)
                    
                    if(!is.na(dominant.veg)){veg.code <- paste(veg.code,dominant.veg,sep="/")} #append dominant understory if there is one
                    #note that the FIA data will produce many tree-vegetation groups that are not listed in the old growth classifications. These will still be able to be assesed for maturity, despite not being assessed for old growth
                    
                  }#end instructions for if vegetation data exists
                  
                  #calculate basal area for plot
                  local.basal.trees <- local.trees[local.trees$DIA >= 5 & local.trees$STATUSCD == 1, ] #isolate to living trees with DBH greater than 5 inches
                  local.basal.trees$basalArea <- (local.basal.trees$DIA^2)*0.005454 #calculate basal area of each tree
                  local.basal.area <- sum(local.basal.trees$basalArea)/condition.area #calculate total basal area per acre of plot
                  
                  #assign old growth type specific to Idaho section of Region 1 based on pages 13-42 of "Old-growth forest types of the Northern Region" by Green, et al. (1992, errata corrected 2011)
                  #note that this is an approximation - actual characterization requires an on-the-ground assessment
                  northern.OG.forest.type <- NA #assign old growth forest type for northern region (this is different from forest.type)
                  if(forest.type %in% c(200:203)){northern.OG.forest.type <- "DF"}
                  if(forest.type %in% c(220:226)){northern.OG.forest.type <- "PP"}
                  if(forest.type %in% c(320,321)){northern.OG.forest.type <- "L"}
                  if(forest.type %in% c(280,281)){northern.OG.forest.type <- "LP"}
                  if(forest.type %in% c(368)){northern.OG.forest.type <- "Y"} #using "miscellaneous western softwoods" to capture pacific yew since no other category exists that doesn't already have a classification in this system
                  if(forest.type %in% 267){northern.OG.forest.type <- "GF"}
                  if(forest.type %in% c(265,266,268)){northern.OG.forest.type <- "SAF"}
                  if(forest.type %in% c(240,241)){northern.OG.forest.type <- "WP"}
                  if(forest.type %in% 301){northern.OG.forest.type <- "WH"}
                  if(forest.type %in% 270){northern.OG.forest.type <- "MAF"}
                  if(forest.type %in% 367){northern.OG.forest.type <- "WBP"}
                  if(forest.type %in% c(304)){northern.OG.forest.type <- "C"}
                  #DAN! WSL??? WHAT IS WSL? SEE WESTERN MONTANA SECTION
                  
                  #determine sub-region of northern section
                  if(locale.state$STUSPS %in% c("ID","WA")){northern.subregion <- "northern Idaho zone"}
                  if(locale.state$STUSPS %in% c("WY","SD","ND")){northern.subregion <- "eastern Montana zone"}
                  cont.divide.test <- st_filter(x = local.plot, #test if plot is east of the continental divide
                                                y = ContDivideEast)
                  if(locale.state$STUSPS %in% "MT" & nrow(cont.divide.test) > 1){northern.subregion <- "eastern Montana zone"}
                  if(locale.state$STUSPS %in% "MT" & nrow(cont.divide.test) < 1){northern.subregion <- "western Montana zone"}
                  
                  
                  if(northern.subregion %in% "northern Idaho zone"){#begin old growth assessment for (northern) region specific to Idaho
                    
                    #assign habitat type group (code) specific to Idaho section of Region 1 (note that Region 1 refers to "PICEA" for many species codes, but this is not a documented code in the plants database. I instead assumed this was in reference to the genus Picea [spruce]) based on Appendix A in "Old-growth forest types of the Northern Region" by Green, et al. (1992, errata corrected 2011)
                    #note that this is an approximation - actual characterization requires an on-the-ground assessment
                    northern.region.codes <- c(NA)
                    if(veg.code %in% c("PIPO","PIPO/AGSP","PIPO/FEID","PIPO/FESC","PIPO/PURT","PIPO/AGSP","PIPO/SYAL","PIPO/BERE","PSME/AGSP","PSME/FEID","PSME/FESC")){northern.region.codes <- c(northern.region.codes,"A")}
                    if(veg.code %in% c("PIPO/PRIV","PIPO/SHCA","PSME/PHMA","PSME/SMST","PSME/CARU","PSME/SYAL","PSME/VAGL","VAGL/ARUV","VAGL/XETE","PSME/VACA","PSME/AGSP","PSME/CARU","PSME/CARU","PSME/AGSP","PSME/ARUV","PSME/PIPO","PSME/CAGE","PSME/SPBE","PSME/ARUV","PSME/JUCO","PSME/ARCO","PSME/SYOR","PIAB/PHMA","PIBR/PHMA","PIEN/PHMA","PIGL/PHMA","PIMA/PHMA","PIPU/PHMA","PIRU/PHMA","PISI/PHMA")){northern.region.codes <- c(northern.region.codes,"B")}
                    if(veg.code %in% c("ABGR/SETR","ABGR/ASCA","ABGR/MEFE","ABGR/ASCA","ABGR/CLUN","ABLA/CLUN","ABLA/ARNU","ABGR/MEFE","ABGR/PHMA")){northern.region.codes <- c(northern.region.codes,"C")}
                    if(veg.code %in% c("ABGR/ASCA","ABGR/TABR","ABGR/CLUN")){northern.region.codes <- c(northern.region.codes,"C1")}
                    if(veg.code %in% c("PSME/LIBO","PSME/SYAL","PSME/CARU","PSME/VAGL","ABGR/LIBO","ABGR/XETE","ABGR/COOC","ABGR/VAGL","ABGR/XETE","ABGR/CLUN","ABGR/VAGL")){northern.region.codes <- c(northern.region.codes,"D")}
                    if(veg.code %in% c("ABGR/PHMA","ABGR/COOC","ABGR/SPBE")){northern.region.codes <- c(northern.region.codes,"E")}
                    if(veg.code %in% c("THPL/OPHO","THPL/ATFI","THPL/ADPE","THPL/ATFI")){northern.region.codes <- c(northern.region.codes,"F")}
                    if(veg.code %in% c("THPL/CLUN","THPL/ARNU","THPL/GYDR","THPL/ASCA","THPL/MEFE","THPL/CLUN","THPL/XETE","TSHE/GYDR","TSHE/ASCA","TSHE/ARNU","TSHE/MEFE","TSHE/ASCA","TSHE/CLUN","TSHE/XETE")){northern.region.codes <- c(northern.region.codes,"G")}
                    if(veg.code %in% c("THPL/CLUN","THPL/TABR","THPL/ASCA")){northern.region.codes <- c(northern.region.codes,"G1")}
                    if(veg.code %in% c("PIAB/EQAR","PIBR/EQAR","PIEN/EQAR","PIGL/EQAR","PIMA/EQAR","PIPU/EQAR","PIRU/EQAR","PISI/EQAR","PIAB/GATR","PIBR/GATR","PIEN/GATR","PIGL/GATR","PIMA/GATR","PIPU/GATR","PIRU/GATR","PISI/GATR","ABLA/OPHO","ABLA/GATR","ABLA/STAM","ABLA/METE","ABLA/LICA","ABLA/CACA","ABLA/LICA","ABLA/GATR","ABLA/VACA","ABLA/LEGL","TSME/STAM","TSME/LUHI","TSME/MEFE")){northern.region.codes <- c(northern.region.codes,"H")}
                    if(veg.code %in% c("PIAB/CLUN","PIBR/CLUN","PIEN/CLUN","PIGL/CLUN","PIMA/CLUN","PIPU/CLUN","PIRU/CLUN","PISI/CLUN","PIAB/LIBO","PIBR/LIBO","PIEN/LIBO","PIGL/LIBO","PIMA/LIBO","PIPU/LIBO","PIRU/LIBO","PISI/LIBO","TSME/MEFE","ABLA/CLUN","ABLA/ARNU","ABLA/VACA","ABLA/XETE","ABLA/MEFE","ABLA/LIBO","ABLA/MEFE","ABLA/COOC","ABLA/LUHI","ABLA/XETE","ABLA/VASC","TSME/CLUN","TSME/MEFE","TSME/XETE")){northern.region.codes <- c(northern.region.codes,"I")}
                    if(veg.code %in% c("PIAB/CLUN","PIBR/CLUN","PIEN/CLUN","PIGL/CLUN","PIMA/CLUN","PIPU/CLUN","PIRU/CLUN","PISI/CLUN","ABLA/VACA","ABLA/XETE","ABLA/VAGL","ABLA/VASC","ABLA/COOC","ABLA/LUHI","TSME/XETE","TSME/LUHI","TSME/VASC","TSME/XETE","ABLA/VAGL","ABLA/THOC","ABLA/VASC","ABLA/CARU","ABLA/CLPS","ABLA/ARCO","ABLA/CAGE","ABLA/PSME")){northern.region.codes <- c(northern.region.codes,"J")}
                    if(veg.code %in% c("ABLA/RIMO","ABLA-PIAL","ABLA/VASC","PIAL/VASC","ABLA/LUHI","ABLA/VASC","TSME/LUHI","ABLA/VASC","ABLA/MEFE","TSME/VASC","TSME/MEFE","PIAL-ABLA","LALY-ABLA","PIAL","PICO/VACA","PICO/XETE","PICO/LIBO","PICO/VASC","PICO/CARU")){northern.region.codes <- c(northern.region.codes,"K")}
                    
                    for(v in 1:length(northern.region.codes)){#begin looping through possible forest groups for (northern) forests in Idaho
                      
                      if(northern.region.codes[v] %in% c("A","B") & northern.OG.forest.type %in% c("PP","DF","L")){#begin sub-assessment for group 1 in (northern) Idaho
                        
                        t.age <- ifelse(raw.stand.age >= 150, 1, 0) #determine if age threshold is met
                        t.density <- ifelse(nrow(local.trees[local.trees$DIA >= 21, ])/condition.area >= 8, 1, 0) #determine if large tree density threshold is met
                        t.basal <- ifelse(local.basal.area >= 40, 1, 0) #determine if basal area threshold is met
                        
                        local.MOG.status <- ifelse(sum(t.age, t.density, t.basal, na.rm=TRUE) == 3, 1, 0) #if all criteria are met, designate plot as old growth (1). Otherwise assign non-old growth (0)
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end sub-assessment for habitat type group 1 (northern) Idaho
                      
                      if(northern.region.codes[v] %in% c("B","C","D","E","G","H","I","J","K") & northern.OG.forest.type %in% c("LP")){#begin sub-assessment for group 2 in (northern) Idaho
                        
                        t.age <- ifelse(raw.stand.age >= 120, 1, 0) #determine if age threshold is met
                        t.density <- ifelse(nrow(local.trees[local.trees$DIA >= 13, ])/condition.area >= 10, 1, 0) #determine if large tree density threshold is met
                        t.basal <- ifelse(local.basal.area >= 60, 1, 0) #determine if basal area threshold is met
                        
                        local.MOG.status <- ifelse(sum(t.age, t.density, t.basal, na.rm=TRUE) == 3, 1, 0) #if all criteria are met, designate plot as old growth (1). Otherwise assign non-old growth (0)
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end sub-assessment for group 2 in (northern) Idaho
                      
                      if(northern.region.codes[v] %in% c("C","C1","G1") & northern.OG.forest.type %in% c("Y")){#begin sub-assessment for habitat type group 3 in (northern) Idaho
                        
                        t.age <- ifelse(raw.stand.age >= 150, 1, 0) #determine if age threshold is met
                        t.density <- ifelse(nrow(local.trees[local.trees$DIA >= 21, ])/condition.area >= 3, 1, 0) #determine if large tree density threshold is met
                        t.basal <- ifelse(local.basal.area >= 80, 1, 0) #determine if basal area threshold is met
                        
                        local.MOG.status <- ifelse(sum(t.age, t.density, t.basal, na.rm=TRUE) == 3, 1, 0) #if all criteria are met, designate plot as old growth (1). Otherwise assign non-old growth (0)
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end sub-assessment for habitat type group 3 in (northern) Idaho
                      
                      if(northern.region.codes[v] %in% c("C","C1","D","E") & northern.OG.forest.type %in% c("DF","GF","SAF","WP","PP")){#begin sub-assessment for habitat type group 4A in (northern) Idaho
                        
                        t.age <- ifelse(raw.stand.age >= 150, 1, 0) #determine if age threshold is met
                        t.density <- ifelse(nrow(local.trees[local.trees$DIA >= 21, ])/condition.area >= 10, 1, 0) #determine if large tree density threshold is met
                        t.basal <- ifelse(local.basal.area >= 80, 1, 0) #determine if basal area threshold is met
                        
                        local.MOG.status <- ifelse(sum(t.age, t.density, t.basal, na.rm=TRUE) == 3, 1, 0) #if all criteria are met, designate plot as old growth (1). Otherwise assign non-old growth (0)
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end sub-assessment for habitat type group 4A in (northern) Idaho
                      
                      if(northern.region.codes[v] %in% c("F","G","G1","H","I") & northern.OG.forest.type %in% c("DF","GF","WH","WP","PP")){#begin sub-assessment for habitat type group 4B in (northern) Idaho
                        
                        t.age <- ifelse(raw.stand.age >= 150, 1, 0) #determine if age threshold is met
                        t.density <- ifelse(nrow(local.trees[local.trees$DIA >= 21, ])/condition.area >= 10, 1, 0) #determine if large tree density threshold is met
                        if(northern.region.codes[v] %in% c("F","G","G1")){t.basal <- ifelse(local.basal.area >= 120, 1, 0)} #determine if basal area threshold is met
                        if(northern.region.codes[v] %in% c("H","I")){t.basal <- ifelse(local.basal.area >= 80, 1, 0)} #determine if basal area threshold is met
                        
                        local.MOG.status <- ifelse(sum(t.age, t.density, t.basal, na.rm=TRUE) == 3, 1, 0) #if all criteria are met, designate plot as old growth (1). Otherwise assign non-old growth (0)
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end sub-assessment for habitat type group 4B in (northern) Idaho
                      
                      if(northern.region.codes[v] %in% c("F","G","G1","H","I") & northern.OG.forest.type %in% c("SAF","MAF")){#begin sub-assessment for habitat type group 5 in (northern) Idaho
                        
                        t.age <- ifelse(raw.stand.age >= 150, 1, 0) #determine if age threshold is met
                        t.density <- ifelse(nrow(local.trees[local.trees$DIA >= 17, ])/condition.area >= 10, 1, 0) #determine if large tree density threshold is met
                        t.basal <- ifelse(local.basal.area >= 80, 1, 0) #determine if basal area threshold is met
                        
                        local.MOG.status <- ifelse(sum(t.age, t.density, t.basal, na.rm=TRUE) == 3, 1, 0) #if all criteria are met, designate plot as old growth (1). Otherwise assign non-old growth (0)
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end sub-assessment for habitat type group 5 in (northern) Idaho
                      
                      if(northern.region.codes[v] %in% c("I","J","K") & northern.OG.forest.type %in% c("WBP")){#begin sub-assessment for habitat type group 6 in (northern) Idaho
                        
                        t.age <- ifelse(raw.stand.age >= 150, 1, 0) #determine if age threshold is met
                        t.density <- ifelse(nrow(local.trees[local.trees$DIA >= 13, ])/condition.area >= 5, 1, 0) #determine if large tree density threshold is met
                        if(northern.region.codes[v] %in% c("I","J")){t.basal <- ifelse(local.basal.area >= 60, 1, 0)} #determine if basal area threshold is met
                        if(northern.region.codes[v] %in% c("K")){t.basal <- ifelse(local.basal.area >= 40, 1, 0)} #determine if basal area threshold is met
                        
                        local.MOG.status <- ifelse(sum(t.age, t.density, t.basal, na.rm=TRUE) == 3, 1, 0) #if all criteria are met, designate plot as old growth (1). Otherwise assign non-old growth (0)
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end sub-assessment for habitat type group 6 in (northern) Idaho
                      
                      if(northern.region.codes[v] %in% c("F","G","G1") & northern.OG.forest.type %in% c("C")){#begin sub-assessment for habitat type group 7 in (northern) Idaho
                        
                        t.age <- ifelse(raw.stand.age >= 150, 1, 0) #determine if age threshold is met
                        
                        t.cedars <- local.trees[local.trees$SPCD %in% c(41,42,43,67,68,81,241,242,978) & local.trees$DIA >= 25, ] #subset to specific species and dbh requirements (specific to this habitat type group)
                        t.specificothers <- local.trees[local.trees$SPCD %in% c(201,202,17,73,263,114,119,122) & local.trees$DIA >= 21, ] #subset to specific species and dbh requirements (specific to this habitat type group)
                        t.others <- local.trees[!(local.trees$SPCD %in% c(41,42,43,67,68,81,241,242,978,201,202,17,73,263,114,119,122)) & local.trees$DIA >= 17, ] #subset to specific species and dbh requirements (specific to this habitat type group)
                        temp.trees <- rbind(t.cedars,t.specificothers,t.others)
                        t.density <- ifelse(nrow(temp.trees)/condition.area >= 10, 1, 0) #determine if large tree density threshold is met
                        
                        t.basal <- ifelse(local.basal.area >= 120, 1, 0) #determine if basal area threshold is met
                        
                        local.MOG.status <- ifelse(sum(t.age, t.density, t.basal, na.rm=TRUE) == 3, 1, 0) #if all criteria are met, designate plot as old growth (1). Otherwise assign non-old growth (0)
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end sub-assessment for habitat type group 7 in (northern) Idaho
                      
                      if(northern.region.codes[v] %in% c("J") & northern.OG.forest.type %in% c("DF","L","SAF","MAF","WP")){#begin sub-assessment for habitat type group 8 in (northern) Idaho
                        
                        t.age <- ifelse(raw.stand.age >= 150, 1, 0) #determine if age threshold is met
                        t.density <- ifelse(nrow(local.trees[local.trees$DIA >= 17, ])/condition.area >= 10, 1, 0) #determine if large tree density threshold is met
                        t.basal <- ifelse(local.basal.area >= 60, 1, 0) #determine if basal area threshold is met
                        
                        local.MOG.status <- ifelse(sum(t.age, t.density, t.basal, na.rm=TRUE) == 3, 1, 0) #if all criteria are met, designate plot as old growth (1). Otherwise assign non-old growth (0)
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end sub-assessment for habitat type group 8 in (northern) Idaho
                      
                      if(northern.region.codes[v] %in% "K" & northern.OG.forest.type %in% c("SAF","MAF")){#begin sub-assessment for habitat type group 9 in (northern) Idaho
                        
                        t.age <- ifelse(raw.stand.age >= 150, 1, 0) #determine if age threshold is met
                        t.density <- ifelse(nrow(local.trees[local.trees$DIA >= 13, ])/condition.area >= 5, 1, 0) #determine if large tree density threshold is met
                        t.basal <- ifelse(local.basal.area >= 40, 1, 0) #determine if basal area threshold is met
                        
                        local.MOG.status <- ifelse(sum(t.age, t.density, t.basal, na.rm=TRUE) == 3, 1, 0) #if all criteria are met, designate plot as old growth (1). Otherwise assign non-old growth (0)
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end sub-assessment for habitat type group 9 in (northern) Idaho
                      
                    }#end looping through possible forest groups for (northern) forests in Idaho
                  }#end old growth assessment for (northern) region specific to Idaho
                  
                  
                  if(northern.subregion %in% "western Montana zone"){#begin old growth assessment for (northern) region specific to west of the Continental Divide
                    
                    #assign habitat type group (code) specific to western Montana section of Region 1 (note that Region 1 refers to "PICEA" for many species codes, but this is not a documented code in the plants database. I instead assumed this was in reference to the genus Picea [spruce]) based on Appendix A in "Old-growth forest types of the Northern Region" by Green, et al. (1992, errata corrected 2011)
                    #note that this is an approximation - actual characterization requires an on-the-ground assessment
                    northern.region.codes <- c(NA)
                    if(veg.code %in% c("PIFL/AGSP","PIFL/FIED","PIFL/FESC","PIFL/JUCO","PIPO/AND","PIPO/AGSP","PIPO/FEID","PIPI/FEID","PIPO/FESC","PIPO/PUTR","PIPO/SYAL","PIPO/PRVI","PIPO/PRIV","PIPO/SHCA","PIPO/PHMA","PSME/AGSP","PSME/FEID","PSME/FESC","PSME/SYAL","PSME/AGSP","PSME/CARU","PSME/ARUV")){northern.region.codes <- c(northern.region.codes,"A")}
                    if(veg.code %in% c("PSME/VACA","PSME/PHMA","PSME/CARU","PSME/SMST","PSME/VAGL","PSME/ARUV","PSME/SYAL","PSME/PIPO","PSME/SPBE","PSME/SYOR","ABGR/SPBE","ABGR/PHMA","ABGR/COOC","ABGR/PHMA")){northern.region.codes <- c(northern.region.codes,"B")}
                    if(veg.code %in% c("PSME/VAGL","PSME/XETE","PSME/LIBO","PSME/CARU","PSME/JUCO","PSME/ARCO","ABGR/XETE","ABGR/COOC","ABGR/VAGL","ABLA/CARU")){northern.region.codes <- c(northern.region.codes,"C")}
                    if(veg.code %in% c("ABGR/VAGL","ABGR/ASCR","ABGR/ASCA","ABGR/MEFE","ABGR/TABR","ABGR/CLUN","ABGR/ARNU","ABGR/XETE","ABGR/PHMA","ABGR/SETR","THPL/CLUN","THPL/ARNU","THPL/MEFE","THPL/XETE","THPL/TABR","THPL/ASCA","THPL/GYDR","TSHE/GYDR","TSHE/CLUN","TSHE/ARNU","TSHE/MEFE","THSE/XETE","TSHE/ASCA")){northern.region.codes <- c(northern.region.codes,"D")}
                    if(veg.code %in% c("PIAB/CLUN","PIBR/CLUN","PIEN/CLUN","PIGL/CLUN","PIMA/CLUN","PIPU/CLUN","PIRU/CLUN","PISI/CLUN","PIOM/CLUN","PIAB/VACA","PIBR/VACA","PIEN/VACA","PIGL/VACA","PIMA/VACA","PIPU/VACA","PIRU/VACA","PISI/VACA","PIOM/VACA","PIAB/SEST","PIBR/SEST","PIEN/SEST","PIGL/SEST","PIMA/SEST","PIPU/SEST","PIRU/SEST","PISI/SEST","PIOM/SEST","PIAB/PSME","PIBR/PSME","PIEN/PSME","PIGL/PSME","PIMA/PSME","PIPU/PSME","PIRU/PSME","PISI/PSME","PIOM/PSME","ABLA/CLUN","ABLA/ARNU","ABLA/VACA","ABLA/XETE","ABLA/MEFE","ABLA/LIBO","ABLA/COOC","ABLA/LUHI","ABLA/VASC","TSME/STAM","TSME/MEFE","TSME/LUHI","TSME/XETE","TSME/CLUN","TSME/MEFE","ABLA/ALSI","ABLA/LUHI","ABLA/MEFE")){northern.region.codes <- c(northern.region.codes,"E")}
                    if(veg.code %in% c("PIAB/EQAR","PIBR/EQAR","PIEN/EQAR","PIGL/EQAR","PIMA/EQAR","PIPU/EQAR","PIRU/EQAR","PISI/EQAR","PIOM/EQAR","PIAB/GATR","PIBR/GATR","PIEN/GATR","PIGL/GATR","PIMA/GATR","PIPU/GATR","PIRU/GATR","PISI/GATR","PIOM/GATR","PIAB/SMST","PIBR/SMST","PIEN/SMST","PIGL/SMST","PIMA/SMST","PIPU/SMST","PIRU/SMST","PISI/SMST","PIOM/SMST","THPL/ALFI","THPL/ADPE","THPL/OPHO","THPL/ADPE","ABLA/OPHO","ABLA/GATR","ABLA/CACA","ABLA/STAM","ABLA/MEFE","ABLA/LICA","ABLA/VACA","ABLA/LEGL","TSME/STAM","TSME/LUHI","TSME/MEFE")){northern.region.codes <- c(northern.region.codes,"F")}
                    if(veg.code %in% c("PSME/LIBO","ABGR/LIBO","ABGR/UBO","PSME/VAGL","ABGR/LIBO")){northern.region.codes <- c(northern.region.codes,"G")}
                    if(veg.code %in% c("PIAB/PHMA","PIBR/PHMA","PIEN/PHMA","PIGL/PHMA","PIMA/PHMA","PIPU/PHMA","PIRU/PHMA","PISI/PHMA","PIOM/PHMA","PIAB/VACA","PIBR/VACA","PIEN/VACA","PIGL/VACA","PIMA/VACA","PIPU/VACA","PIRU/VACA","PISI/VACA","PIOM/VACA","ABLA/VACA","ABLA/LIBO","ABLA/VASC","ABLA/XETE","ABLA/VAGL","ABLA/COOC","ABLA/LUHI","TSME/XETE","TSME/LUHI","TSME/VAGL","ABLA/VAGL","ABLA/VASC","ABLA/THOC","ABLA/ARCO","ABLA/CAGE","ABLA/GAGE","ABLA/PSME","ABLA/RIMO","ABAL/RIMO","PICO/PUTR","PICO/VACA","PICO/XETE","PICO/LIBO","PICO/VASC","PICO/CARU")){northern.region.codes <- c(northern.region.codes,"H")}
                    if(veg.code %in% c("ABLA/VASC","ABLA/PIAL","ABLA/LUHI","TSME/LUHI","TSME/VASC","TSME/MEFE")){northern.region.codes <- c(northern.region.codes,"I")}
                    if(veg.code %in% c("PIAL/ABLA","PIAL","LALY/ABLA")){northern.region.codes <- c(northern.region.codes,"J")}
                    
                    
                    for(v in 1:length(northern.region.codes)){#begin looping through possible forest groups for (northern) forests in western Montana
                      
                      if(northern.region.codes[v] %in% c("A","B") & northern.OG.forest.type %in% c("PP","DF","L","GF","LP")){#begin sub-assessment for group 1 in (northern) western Montana
                        
                        t.age <- ifelse(raw.stand.age >= 170, 1, 0) #determine if age threshold is met
                        t.density <- ifelse(nrow(local.trees[local.trees$DIA >= 21, ])/condition.area >= 8, 1, 0) #determine if large tree density threshold is met
                        t.basal <- ifelse(local.basal.area >= 60, 1, 0) #determine if basal area threshold is met
                        
                        local.MOG.status <- ifelse(sum(t.age, t.density, t.basal, na.rm=TRUE) == 3, 1, 0) #if all criteria are met, designate plot as old growth (1). Otherwise assign non-old growth (0)
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end sub-assessment for habitat type group 1 (northern) western Montana
                      
                      if(northern.region.codes[v] %in% c("C") & northern.OG.forest.type %in% c("DF","L","PP","SAF","GF")){#begin sub-assessment for group 2 in (northern) western Montana
                        
                        t.age <- ifelse(raw.stand.age >= 170, 1, 0) #determine if age threshold is met
                        t.density <- ifelse(nrow(local.trees[local.trees$DIA >= 21, ])/condition.area >= 8, 1, 0) #determine if large tree density threshold is met
                        t.basal <- ifelse(local.basal.area >= 80, 1, 0) #determine if basal area threshold is met
                        
                        local.MOG.status <- ifelse(sum(t.age, t.density, t.basal, na.rm=TRUE) == 3, 1, 0) #if all criteria are met, designate plot as old growth (1). Otherwise assign non-old growth (0)
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end sub-assessment for habitat type group 2 (northern) western Montana
                      
                      if(northern.region.codes[v] %in% c("C","D","E","F","G","H") & northern.OG.forest.type %in% c("LP")){#begin sub-assessment for group 3 in (northern) western Montana
                        
                        t.age <- ifelse(raw.stand.age >= 140, 1, 0) #determine if age threshold is met
                        t.density <- ifelse(nrow(local.trees[local.trees$DIA >= 13, ])/condition.area >= 10, 1, 0) #determine if large tree density threshold is met
                        t.basal <- ifelse(local.basal.area >= 80, 1, 0) #determine if basal area threshold is met
                        if(northern.region.codes[v] %in% "E"){t.basal <- ifelse(local.basal.area >= 60, 1, 0)} #determine if basal area threshold is met
                        if(northern.region.codes[v] %in% c("C","H")){t.basal <- ifelse(local.basal.area >= 70, 1, 0)} #determine if basal area threshold is met
                        
                        local.MOG.status <- ifelse(sum(t.age, t.density, t.basal, na.rm=TRUE) == 3, 1, 0) #if all criteria are met, designate plot as old growth (1). Otherwise assign non-old growth (0)
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end sub-assessment for habitat type group 3 (northern) western Montana
                      
                      if(northern.region.codes[v] %in% c("D","E","F") & northern.OG.forest.type %in% c("SAF","DF","GF","C","L","MAF","PP","WP","WH","WSL")){#begin sub-assessment for group 4 in (northern) western Montana
                        
                        t.age <- ifelse(raw.stand.age >= 180, 1, 0) #determine if age threshold is met
                        t.density <- ifelse(nrow(local.trees[local.trees$DIA >= 21, ])/condition.area >= 10, 1, 0) #determine if large tree density threshold is met
                        t.basal <- ifelse(local.basal.area >= 80, 1, 0) #determine if basal area threshold is met
                        
                        local.MOG.status <- ifelse(sum(t.age, t.density, t.basal, na.rm=TRUE) == 3, 1, 0) #if all criteria are met, designate plot as old growth (1). Otherwise assign non-old growth (0)
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end sub-assessment for habitat type group 4 (northern) western Montana
                      
                      if(northern.region.codes[v] %in% c("G","H") & northern.OG.forest.type %in% c("SAF","DF","GF","L","MAF","PP","WP","WSL")){#begin sub-assessment for group 5 in (northern) western Montana
                        
                        t.age <- ifelse(raw.stand.age >= 180, 1, 0) #determine if age threshold is met
                        t.density <- ifelse(nrow(local.trees[local.trees$DIA >= 17, ])/condition.area >= 10, 1, 0) #determine if large tree density threshold is met
                        t.basal <- ifelse(local.basal.area >= 80, 1, 0) #determine if basal area threshold is met
                        if(northern.region.codes[v] %in% "H" & northern.OG.forest.type %in% "SAF"){t.basal <- ifelse(local.basal.area >= 70, 1, 0)}
                        
                        local.MOG.status <- ifelse(sum(t.age, t.density, t.basal, na.rm=TRUE) == 3, 1, 0) #if all criteria are met, designate plot as old growth (1). Otherwise assign non-old growth (0)
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end sub-assessment for habitat type group 5 (northern) western Montana
                      
                      if(northern.region.codes[v] %in% c("I") & northern.OG.forest.type %in% c("SAF","WSL","DF","L")){#begin sub-assessment for group 6 in (northern) western Montana
                        
                        t.age <- ifelse(raw.stand.age >= 180, 1, 0) #determine if age threshold is met
                        t.density <- ifelse(nrow(local.trees[local.trees$DIA >= 13, ])/condition.area >= 10, 1, 0) #determine if large tree density threshold is met
                        t.basal <- ifelse(local.basal.area >= 60, 1, 0) #determine if basal area threshold is met
                        
                        local.MOG.status <- ifelse(sum(t.age, t.density, t.basal, na.rm=TRUE) == 3, 1, 0) #if all criteria are met, designate plot as old growth (1). Otherwise assign non-old growth (0)
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end sub-assessment for habitat type group 6 (northern) western Montana
                      
                      if(northern.region.codes[v] %in% c("I") & northern.OG.forest.type %in% c("LP")){#begin sub-assessment for group 7 in (northern) western Montana
                        
                        t.age <- ifelse(raw.stand.age >= 140, 1, 0) #determine if age threshold is met
                        t.density <- ifelse(nrow(local.trees[local.trees$DIA >= 9, ])/condition.area >= 30, 1, 0) #determine if large tree density threshold is met
                        t.basal <- ifelse(local.basal.area >= 70, 1, 0) #determine if basal area threshold is met
                        
                        local.MOG.status <- ifelse(sum(t.age, t.density, t.basal, na.rm=TRUE) == 3, 1, 0) #if all criteria are met, designate plot as old growth (1). Otherwise assign non-old growth (0)
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end sub-assessment for habitat type group 7 (northern) western Montana
                      
                      if(northern.region.codes[v] %in% c("J") & northern.OG.forest.type %in% c("SAF","WSL")){#begin sub-assessment for group 8 in (northern) western Montana
                        
                        t.age <- ifelse(raw.stand.age >= 180, 1, 0) #determine if age threshold is met
                        t.density <- ifelse(nrow(local.trees[local.trees$DIA >= 13, ])/condition.area >= 20, 1, 0) #determine if large tree density threshold is met
                        t.basal <- ifelse(local.basal.area >= 80, 1, 0) #determine if basal area threshold is met
                        
                        local.MOG.status <- ifelse(sum(t.age, t.density, t.basal, na.rm=TRUE) == 3, 1, 0) #if all criteria are met, designate plot as old growth (1). Otherwise assign non-old growth (0)
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end sub-assessment for habitat type group 8 (northern) western Montana
                      
                    }#end looping through possible forest groups for (northern) forests in western Montana
                  }#end old growth assessment for (northern) region specific to west of the Continental Divide
                  
                  if(northern.subregion %in% "eastern Montana zone"){#begin old growth assessment for (northern) region specific to west of the Continental Divide
                    
                    #assign habitat type group (code) specific to western Montana section of Region 1 (note that Region 1 refers to "PICEA" for many species codes, but this is not a documented code in the plants database. I instead assumed this was in reference to the genus Picea [spruce]) based on Appendix A in "Old-growth forest types of the Northern Region" by Green, et al. (1992, errata corrected 2011)
                    #note that this is an approximation - actual characterization requires an on-the-ground assessment
                    northern.region.codes <- c(NA)
                    if(veg.code %in% c("PIFL","PIFL/AGSP","PIFL/FEID","PIFL/FESC","PIFL/JUCO","PIPO","PIPO/JUSC","PIPO/MARE","PIPO/AND","PIPO/AGSP","PIPO/FEID","PIPO/FESC","PIPO/SYOC","JUSC/PSSP","JUSC/PSSPS","PIPO/PUTR","PIPO/AGSP","PIPO/FEID","PIPO/FESC","PIPO/SYOC","JUSC/PSSP","JUSC/PSSPS","PIPO/PUTR","PIPO/AGSP","PIPO/FEID","PIPO/SYAL","PIPO/PRVI","PIPO/SHCA","PIPO/PHMA","PSME/AGSP","PSME/FEID","PSME/FESC","PSME/SYAL","PSME/AGSP","PSME/CARU","PSME/AGSP","PSME/ARUV","PSME/JUCO","PSME/ARCO","PSME/SYOR","PSME/JUSC")){northern.region.codes <- c(northern.region.codes,"A")}
                    if(veg.code %in% c("PSME/PHMA","PSME/SMST","PSME/CARU","PSME/ARUV","PSME/PIPO","PSME/CAGE","PSME/SPBE","PICO/PUTR")){northern.region.codes <- c(northern.region.codes,"B")}
                    if(veg.code %in% c("PSME/MARE","PSME/VACA","PSME/PHMA","PSME/VAGL","PSME/ARUV","PSME/PHMA","PSME/CARU","PSME/SYAL","PSME/CARU","PIAB/PHMA","PIBR/PHMA","PIEN/PHMA","PIGL/PHMA","PIMA/PHMA","PIPU/PHMA","PIRU/PHMA","PISI/PHMA","PIOM/PHMA")){northern.region.codes <- c(northern.region.codes,"C")}
                    if(veg.code %in% c("PSME/VAGL","PSME/XETE","PSME/LIBO","PSME/SYAL","PSME/CARU","PSME/VAGL","PIAB/LIBO","PIBR/LIBO","PIEN/LIBO","PIGL/LIBO","PIMA/LIBO","PIPU/LIBO","PIRU/LIBO","PISI/LIBO","PIOM/LIBO","PIAB/SMST","PIBR/SMST","PIEN/SMST","PIGL/SMST","PIMA/SMST","PIPU/SMST","PIRU/SMST","PISI/SMST","PIOM/SMST","ABLA/LIBO","ABLA/XETE","ABLA/VASC","PICO/LIBO")){northern.region.codes <- c(northern.region.codes,"D")}
                    if(veg.code %in% c("PIAB/EQAR","PIBR/EQAR","PIEN/EQAR","PIGL/EQAR","PIMA/EQAR","PIPU/EQAR","PIRU/EQAR","PISI/EQAR","PIOM/EQAR","PIAB/CLUN","PIBR/CLUN","PIEN/CLUN","PIGL/CLUN","PIMA/CLUN","PIPU/CLUN","PIRU/CLUN","PISI/CLUN","PIOM/CLUN","PIAB/VACA","PIBR/VACA","PIEN/VACA","PIGL/VACA","PIMA/VACA","PIPU/VACA","PIRU/VACA","PISI/VACA","PIOM/VACA","PIAB/GATR","PIBR/GATR","PIEN/GATR","PIGL/GATR","PIMA/GATR","PIPU/GATR","PIRU/GATR","PISI/GATR","PIOM/GATR","ABLA/CLUN","ABLA/ARNU","ABLA/VACA","ABLA/XETE","ABLA/MEFE","ABLA/GATR","ALBL/GATR","ALBL/VASC","ABLA/CACA","ABLA/GATR","ABLA/VACA","ABLA/LEGL")){northern.region.codes <- c(northern.region.codes,"E")}
                    if(veg.code %in% c("PIAB/VACA","PIBR/VACA","PIEN/VACA","PIGL/VACA","PIMA/VACA","PIPU/VACA","PIRU/VACA","PISI/VACA","PIOM/VACA","ABLA/SYAL","ABLA/VACA","ABLA/XETE","ABLA/VAGL","ABLA/VASC","ABLA/XETE","TSME/XETE","ABLA/VAGL","ABLA/VASC","ABLA/CARU","ABLA/THOC","PICO/VACA","PICO/VASC","PICO/CARU","PICO/JUCO")){northern.region.codes <- c(northern.region.codes,"F")}
                    if(veg.code %in% c("ABLA/MEFE","ABLA/ALSI")){northern.region.codes <- c(northern.region.codes,"G")}
                    if(veg.code %in% c("PIAB/SEST","PIBR/SEST","PIEN/SEST","PIGL/SEST","PIMA/SEST","PIPU/SEST","PIRU/SEST","PISI/SEST","PIOM/SEST","PICA/JUCO","ABLA/JUCO","ABLA/CARU","ABLA/CLPS","PIAB/PSME","PIBR/PSME","PIEN/PSME","PIGL/PSME","PIMA/PSME","PIPU/PSME","PIRU/PSME","PISI/PSME","PIOM/PSME","ABLA/ARCO","ABLA/CAGE","ABLA/PSME")){northern.region.codes <- c(northern.region.codes,"H")}
                    if(veg.code %in% c("TSME/MEFE","ABLA/RIMO","ABLA/VASC","ABLA'LUHI","ABLA/MEFE","TSME/LUHI","TSME/VASC","TSME/XETE")){northern.region.codes <- c(northern.region.codes,"I")}
                    if(veg.code %in% c("PIAL-ABLA","LALY-ABLA","PIAL")){northern.region.codes <- c(northern.region.codes,"J")}
                    
                    for(v in 1:length(northern.region.codes)){#begin looping through possible forest groups for (northern) forests in eastern Montana
                      
                      if(northern.region.codes[v] %in% c("A") & northern.OG.forest.type %in% c("DF")){#begin sub-assessment for group 1 in (northern) eastern Montana
                        
                        t.age <- ifelse(raw.stand.age >= 200, 1, 0) #determine if age threshold is met
                        t.density <- ifelse(nrow(local.trees[local.trees$DIA >= 17, ])/condition.area >= 4, 1, 0) #determine if large tree density threshold is met
                        t.basal <- ifelse(local.basal.area >= 60, 1, 0) #determine if basal area threshold is met
                        
                        local.MOG.status <- ifelse(sum(t.age, t.density, t.basal, na.rm=TRUE) == 3, 1, 0) #if all criteria are met, designate plot as old growth (1). Otherwise assign non-old growth (0)
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end sub-assessment for habitat type group 1 (northern) eastern Montana
                      
                      if(northern.region.codes[v] %in% c("B","C","D","E","F","H") & northern.OG.forest.type %in% c("DF")){#begin sub-assessment for group 2 in (northern) eastern Montana
                        
                        t.age <- ifelse(raw.stand.age >= 200, 1, 0) #determine if age threshold is met
                        t.density <- ifelse(nrow(local.trees[local.trees$DIA >= 19, ])/condition.area >= 5, 1, 0) #determine if large tree density threshold is met
                        t.basal <- ifelse(local.basal.area >= 60, 1, 0) #determine if basal area threshold is met
                        
                        local.MOG.status <- ifelse(sum(t.age, t.density, t.basal, na.rm=TRUE) == 3, 1, 0) #if all criteria are met, designate plot as old growth (1). Otherwise assign non-old growth (0)
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end sub-assessment for habitat type group 2 (northern) eastern Montana
                      
                      if(northern.region.codes[v] %in% c("G") & northern.OG.forest.type %in% c("DF")){#begin sub-assessment for group 3 in (northern) eastern Montana
                        
                        t.age <- ifelse(raw.stand.age >= 180, 1, 0) #determine if age threshold is met
                        t.density <- ifelse(nrow(local.trees[local.trees$DIA >= 17, ])/condition.area >= 10, 1, 0) #determine if large tree density threshold is met
                        t.basal <- ifelse(local.basal.area >= 80, 1, 0) #determine if basal area threshold is met
                        
                        local.MOG.status <- ifelse(sum(t.age, t.density, t.basal, na.rm=TRUE) == 3, 1, 0) #if all criteria are met, designate plot as old growth (1). Otherwise assign non-old growth (0)
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end sub-assessment for habitat type group 3 (northern) eastern Montana
                      
                      if(northern.region.codes[v] %in% c("A","B","C","K") & northern.OG.forest.type %in% c("PP")){#begin sub-assessment for group 4 in (northern) eastern Montana
                        
                        t.age <- ifelse(raw.stand.age >= 180, 1, 0) #determine if age threshold is met
                        t.density <- ifelse(nrow(local.trees[local.trees$DIA >= 17, ])/condition.area >= 4, 1, 0) #determine if large tree density threshold is met
                        t.basal <- ifelse(local.basal.area >= 40, 1, 0) #determine if basal area threshold is met
                        
                        local.MOG.status <- ifelse(sum(t.age, t.density, t.basal, na.rm=TRUE) == 3, 1, 0) #if all criteria are met, designate plot as old growth (1). Otherwise assign non-old growth (0)
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end sub-assessment for habitat type group 4 (northern) eastern Montana
                      
                      if(northern.region.codes[v] %in% c("A","B") & northern.OG.forest.type %in% c("PF")){#begin sub-assessment for group 5 in (northern) eastern Montana
                        
                        t.age <- ifelse(raw.stand.age >= 120, 1, 0) #determine if age threshold is met
                        t.density <- ifelse(nrow(local.trees[local.trees$DIA >= 9, ])/condition.area >= 6, 1, 0) #determine if large tree density threshold is met
                        t.basal <- ifelse(local.basal.area >= 50, 1, 0) #determine if basal area threshold is met
                        
                        local.MOG.status <- ifelse(sum(t.age, t.density, t.basal, na.rm=TRUE) == 3, 1, 0) #if all criteria are met, designate plot as old growth (1). Otherwise assign non-old growth (0)
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end sub-assessment for habitat type group 5 (northern) eastern Montana
                      
                      if(northern.region.codes[v] %in% c("A","B","C","D","E","F","G","H","I") & northern.OG.forest.type %in% c("LP")){#begin sub-assessment for group 6 in (northern) eastern Montana
                        
                        t.age <- ifelse(raw.stand.age >= 150, 1, 0) #determine if age threshold is met
                        t.density <- ifelse(nrow(local.trees[local.trees$DIA >= 10, ])/condition.area >= 12, 1, 0) #determine if large tree density threshold is met
                        t.basal <- ifelse(local.basal.area >= 50, 1, 0) #determine if basal area threshold is met
                        
                        local.MOG.status <- ifelse(sum(t.age, t.density, t.basal, na.rm=TRUE) == 3, 1, 0) #if all criteria are met, designate plot as old growth (1). Otherwise assign non-old growth (0)
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end sub-assessment for habitat type group 6 (northern) eastern Montana
                      
                      if(northern.region.codes[v] %in% c("C") & northern.OG.forest.type %in% c("SAF")){#begin sub-assessment for group 7 in (northern) eastern Montana
                        
                        t.age <- ifelse(raw.stand.age >= 160, 1, 0) #determine if age threshold is met
                        t.density <- ifelse(nrow(local.trees[local.trees$DIA >= 17, ])/condition.area >= 12, 1, 0) #determine if large tree density threshold is met
                        t.basal <- ifelse(local.basal.area >= 80, 1, 0) #determine if basal area threshold is met
                        
                        local.MOG.status <- ifelse(sum(t.age, t.density, t.basal, na.rm=TRUE) == 3, 1, 0) #if all criteria are met, designate plot as old growth (1). Otherwise assign non-old growth (0)
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end sub-assessment for habitat type group 7 (northern) eastern Montana
                      
                      if(northern.region.codes[v] %in% c("D","E") & northern.OG.forest.type %in% c("SAF")){#begin sub-assessment for group 8 in (northern) eastern Montana
                        
                        t.age <- ifelse(raw.stand.age >= 160, 1, 0) #determine if age threshold is met
                        t.density <- ifelse(nrow(local.trees[local.trees$DIA >= 17, ])/condition.area >= 7, 1, 0) #determine if large tree density threshold is met
                        t.basal <- ifelse(local.basal.area >= 80, 1, 0) #determine if basal area threshold is met
                        
                        local.MOG.status <- ifelse(sum(t.age, t.density, t.basal, na.rm=TRUE) == 3, 1, 0) #if all criteria are met, designate plot as old growth (1). Otherwise assign non-old growth (0)
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end sub-assessment for habitat type group 8 (northern) eastern Montana
                      
                      if(northern.region.codes[v] %in% c("F","G","H","I") & northern.OG.forest.type %in% c("SAF")){#begin sub-assessment for group 9 in (northern) eastern Montana
                        
                        t.age <- ifelse(raw.stand.age >= 160, 1, 0) #determine if age threshold is met
                        t.density <- ifelse(nrow(local.trees[local.trees$DIA >= 13, ])/condition.area >= 10, 1, 0) #determine if large tree density threshold is met
                        t.basal <- ifelse(local.basal.area >= 60, 1, 0) #determine if basal area threshold is met
                        
                        local.MOG.status <- ifelse(sum(t.age, t.density, t.basal, na.rm=TRUE) == 3, 1, 0) #if all criteria are met, designate plot as old growth (1). Otherwise assign non-old growth (0)
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end sub-assessment for habitat type group 9 (northern) eastern Montana
                      
                      if(northern.region.codes[v] %in% c("J") & northern.OG.forest.type %in% c("SAF")){#begin sub-assessment for group 10 in (northern) eastern Montana
                        
                        t.age <- ifelse(raw.stand.age >= 135, 1, 0) #determine if age threshold is met
                        t.density <- ifelse(nrow(local.trees[local.trees$DIA >= 13, ])/condition.area >= 8, 1, 0) #determine if large tree density threshold is met
                        t.basal <- ifelse(local.basal.area >= 40, 1, 0) #determine if basal area threshold is met
                        
                        local.MOG.status <- ifelse(sum(t.age, t.density, t.basal, na.rm=TRUE) == 3, 1, 0) #if all criteria are met, designate plot as old growth (1). Otherwise assign non-old growth (0)
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end sub-assessment for habitat type group 10 (northern) eastern Montana
                      
                      if(northern.region.codes[v] %in% c("D","E","F","G","H","I") & northern.OG.forest.type %in% c("WBP")){#begin sub-assessment for group 11 in (northern) eastern Montana
                        
                        t.age <- ifelse(raw.stand.age >= 150, 1, 0) #determine if age threshold is met
                        t.density <- ifelse(nrow(local.trees[local.trees$DIA >= 13, ])/condition.area >= 11, 1, 0) #determine if large tree density threshold is met
                        t.basal <- ifelse(local.basal.area >= 60, 1, 0) #determine if basal area threshold is met
                        
                        local.MOG.status <- ifelse(sum(t.age, t.density, t.basal, na.rm=TRUE) == 3, 1, 0) #if all criteria are met, designate plot as old growth (1). Otherwise assign non-old growth (0)
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end sub-assessment for habitat type group 11 (northern) eastern Montana
                      
                      if(northern.region.codes[v] %in% c("J") & northern.OG.forest.type %in% c("WBP")){#begin sub-assessment for group 12 in (northern) eastern Montana
                        
                        t.age <- ifelse(raw.stand.age >= 135, 1, 0) #determine if age threshold is met
                        t.density <- ifelse(nrow(local.trees[local.trees$DIA >= 13, ])/condition.area >= 7, 1, 0) #determine if large tree density threshold is met
                        t.basal <- ifelse(local.basal.area >= 40, 1, 0) #determine if basal area threshold is met
                        
                        local.MOG.status <- ifelse(sum(t.age, t.density, t.basal, na.rm=TRUE) == 3, 1, 0) #if all criteria are met, designate plot as old growth (1). Otherwise assign non-old growth (0)
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end sub-assessment for habitat type group 12 (northern) eastern Montana
                      
                    }#end looping through possible forest groups for (northern) forests in eastern Montana
                  }#end old growth assessment for (northern) region specific to east of the Continental Divide
                  
                  #if the plot is not old growth, test to see if it is at least mature based on Table 19 of MOG document
                  #note that, unlike other regions, we assess for maturity automatically in the northern region since each plot can fall into several different categories, so no single value for the plot can be relied on
                  if(northern.OG.forest.type %in% c("DF") | forest.type %in% c(200:203)){#begin mature assessment for (northern) douglas fir
                    
                    t.ddiscore <- ifelse(ddiscore >= 32.6, 0.34, 0)
                    t.badom <- ifelse(badom >= 82.5, 0.33, 0)
                    t.QMDdom <- ifelse(QMDdom >= 10.3, 0.33, 0)
                    
                    local.MOG.status <- sum(t.ddiscore,t.badom,t.QMDdom,na.rm=TRUE)
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end mature assessment for (northern) douglas fir
                  
                  if(northern.OG.forest.type %in% c("SAF","GF","WP") | forest.type %in% c(260:271)){#begin mature assessment for (northern) fir/spruce/mountain hemlock
                    
                    t.ddiscore <- ifelse(ddiscore >= 24, 0.44, 0)
                    t.HTsd <- ifelse(HTsd >= 49.6, 0.3, 0)
                    t.HTquart <- ifelse(HTquart >= 39.2, 0.26, 0)
                    
                    local.MOG.status <- sum(t.ddiscore,t.HTsd,t.HTquart,na.rm=TRUE)
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end mature assessment for (northern) fir/spruce/mountain hemlock
                  
                  if(forest.type %in% c(910:912,700:722,900:905,500:520,970:976)){#begin mature assessment for (northern) hardwoods
                    
                    t.ddiscore <- ifelse(ddiscore >= 23.9, 0.31, 0)
                    t.badom <- ifelse(badom >= 62, 0.28, 0)
                    t.HTquart <- ifelse(HTquart >= 38.4, 0.26, 0)
                    t.HTsd <- ifelse(HTsd >= 28, 0.15, 0)
                    
                    local.MOG.status <- sum(t.ddiscore,t.badom,t.HTsd,t.HTquart,na.rm=TRUE)
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end mature assessment for (northern) hardwoods
                  
                  if(forest.type %in% c(300:305)){#begin mature assessment for (northern) hemlock/sitka spruce
                    
                    t.ddiscore <- ifelse(ddiscore >= 45, 0.38, 0)
                    t.HTsd <- ifelse(HTsd >= 74.4, 0.28, 0)
                    t.HTquart <- ifelse(HTquart >= 69.2, 0.21, 0)
                    t.tpadom <- ifelse(tpadom <= 70, 0.13, 0)
                    
                    local.MOG.status <- sum(t.ddiscore,t.HTsd,t.HTquart,t.tpadom,na.rm=TRUE)
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end mature assessment for (northern) hemlock/sitka spruce
                  
                  if(northern.OG.forest.type %in% c("LP") | forest.type %in% c(280,281)){#begin mature assessment for (northern) lodgepole pine
                    
                    t.HTquart <- ifelse(HTquart >= 25, 0.28, 0)
                    t.ddiscore <- ifelse(ddiscore >= 14.6, 0.26, 0)
                    t.badom <- ifelse(badom >= 43.6, 0.26, 0)
                    t.HTsd <- ifelse(HTsd >= 24, 0.19, 0)
                    
                    local.MOG.status <- sum(t.HTquart,t.ddiscore,t.badom,t.HTsd,na.rm=TRUE)
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end mature assessment for (northern) lodgepole pine
                  
                  if(forest.type %in% c(360:369,180:185)){#begin mature assessment for (northern) pinyon juniper/western softwoods
                    
                    t.ddiscore <- ifelse(ddiscore >- 24, 0.3, 0)
                    t.HTquart <- ifelse(HTquart >= 28.6, 0.25, 0)
                    t.QMDdom <- ifelse(QMDdom >= 7, 0.25, 0)
                    t.HTsd <- ifelse(HTsd >= 29.4, 0.2, 0)
                    
                    local.MOG.status <- sum(t.ddiscore,t.HTquart,t.QMDdom,t.HTsd,na.rm=TRUE)
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end mature assessment for (northern) pinyon juniper/western softwoods
                  
                  if(northern.OG.forest.type %in% c("PP") | forest.type %in% c(220:226)){#begin mature assessment for (northern) ponderosa pine
                    
                    t.ddiscore <- ifelse(ddiscore >= 31.5, 0.36, 0)
                    t.QMDdom <- ifelse(QMDdom >= 13, 0.34, 0)
                    t.HTsd <- ifelse(HTsd >= 40.7, 0.3, 0)
                    
                    local.MOG.status <- sum(t.ddiscore,t.QMDdom,t.HTsd,na.rm=TRUE)
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end mature assessment for (northern) ponderosa pine
                  
                  if(northern.OG.forest.type %in% c("L")){#begin mature assessment for (northern) ponderosa pine
                    
                    t.QMDdom <- ifelse(QMDdom >= 15.8, 0.31, 0)
                    t.ddiscore <- ifelse(ddiscore >= 53, 0.31, 0)
                    t.HTsd <- ifelse(HTsd >= 80.9, 0.21, 0)
                    t.tpadom <- ifelse(tpadom <= 69, 0.16, 0)
                    
                    local.MOG.status <- sum(t.QMDdom,t.ddiscore,t.HTsd,t.tpadom,na.rm=TRUE)
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end mature assessment for (northern) ponderosa pine
                  
                }#end instructions for northern region
                
                if(region %in% "southwest"){#begin instructions for southwest region
                  
                  #define ecological response unit per Table 10 in the MOG document
                  ERU <- "other" #assign as NA to catch errors downstream
                  if(local.condition$HABTYPCD1 %in% c(415,435,604,1100,3060,3080,3090,3110,3111,3112,3200,3201,3202,3203,3231,3240,3300,3301,3310,3320,3350,3370,3999,4060,4061,4062,4151,4152,4300,4310,4320,4330,4340,4350,4351,4360,4999,26005,240300)){ERU <- "spruce-fir forest"}
                  if(local.condition$HABTYPCD1 %in% c(1010,1011,1012,1020,1030,1070,1080,1081,1110,1111,1120,1150,1160,1231,1999,6010,6060,6070,6071,6080,6130,12320,12333)){ERU <- "mixed conifer with aspen"}
                  if(local.condition$HABTYPCD1 %in% c(238040,238310)){ERU <- "bristlecone pine"}
                  if(local.condition$HABTYPCD1 %in% c(1021,1022,1040,1041,1042,1050,1051,1052,1053,1054,1060,1090,1140,1141,1203,1213,1239,1241,6090,11130,12140,12141,12142,12143,12330,12331,12332,12340,12341,12350,12360,12361,12362,12380,12420,12430,12999,238300)){ERU <- "mixed conifer - frequent fire"}
                  if(local.condition$HABTYPCD1 %in% c(11030,11031,11032,11033,11035,11090,11091,11092,11093,11210,11211,11212,11213,11214,11215,11216,11320,11330,11340,11341,11350,11380,11390,11391,11392,11400,11460,11500,11999)){ERU <- "ponderosa pine forest"}
                  if(local.condition$HABTYPCD1 %in% c(11034,11220,11360,11361,11370,11410,11411,11420,11430,11440,32010,32030,32999,33010,33020,33030)){ERU <- "ponderosa pine - evergreen oak"}
                  if(local.condition$HABTYPCD1 %in% c(3102,204400,230030,230040,230041,230042,230999,231010,232070,233010,233030,233040,233041,233042,233050)){ERU <- "pinyon juniper evergreen shrub"}
                  if(local.condition$HABTYPCD1 %in% c(202500,202500,204320,204330,204500,232020,232330,233330)){ERU <- "pinyon juniper (persistent)"}
                  if(local.condition$HABTYPCD1 %in% c(20406,20410,20411,20431,23204,204021,204022,204023,204024,204300,204350,204370,204999,231020,232030,232999,233020,233021,233022,233999,9000042)){ERU <- "pinyon juniper sagebrush"}
                  if(local.condition$HABTYPCD1 %in% c(20404,204050,204321,2040303)){ERU <- "pinyon juniper deciduous shrub"}
                  if(local.condition$HABTYPCD1 %in% c(20406,20410,20411,20431,23204,204021,204022,204023,204024,204300,204350,204370,204999,231020,232030,232999,233020,233021,233022,233999,9000042)){ERU <- "pinyon juniper grass"}
                  if(local.condition$HABTYPCD1 %in% c(20140,201010,201011,201020,201040,201331,201332,201333,201340,201350,201400,201410,201999,202320,202321,202330,202331,202999,231021,231030,231040,231050,231999,9000043)){ERU <- "juniper grass"}
                  if(local.condition$HABTYPCD1 %in% c(3101,204360,232050,232060,630010,630030,630040,630043,630050,2040301,2040302)){ERU <- "madrean pinyon-oak"}
                  if(local.condition$HABTYPCD1 %in% c(31999,610010,610020,620010,620020,620021,620030,620999,630020,630041,630042,632999,650010,650999)){ERU <- "madrean encinal woodland"}
                  if(local.condition$HABTYPCD1 %in% 640999){ERU <- "gambel oak shrubland"}
                  if(local.condition$HABTYPCD1 %in% c(201420, 201430, 210999)){ERU <- "semi-desert grassland"}
                  if(local.condition$HABTYPCD1 %in% 11470){ERU <- "ponderosa pine/willow"}
                  if(local.condition$HABTYPCD1 %in% c(1130, 620040)){ERU <- "arizona walnut"}
                  if(local.condition$HABTYPCD1 %in% 104){ERU <- "rio grande cottonwood/shrub"}
                  if(local.condition$HABTYPCD1 %in% 103){ERU <- "narrowleaf cottonwood-spruce, narrowleaf cottonwood/shrub"}
                  if(local.condition$HABTYPCD1 %in% 3){ERU <- "upper montane conifer/willow"}
                  
                  #calculate Zeide's stand density index 
                  #note that although the official old growth definition for this region requires this calculation, there is little guidance on how to measure it. Given SW-specific citations do not mention this metric. Zeide 1983 offers some guidance, but fails to explain the constants employed in the equation. As a result, I do the best I can here.
                  Zeide.b <- 1.6064 #constant specified in Zeide 1983
                  calc.trees <- local.trees[!is.na(local.trees$DIA), ] #remove trees without DBH data
                  Dr_all <- ((1/nrow(calc.trees))*sum(calc.trees$DIA^Zeide.b))^(1/Zeide.b) #calculate Reineke's diameter for all trees in condition using Equation 3 of Zeide 1983, aseqation 5 requires an unknown coefficient of variation but achieves the same outcome
                  N_all <- Dr_all^-Zeide.b #use modified equation 1 from Zeide 1983 to arrive at relative density (a is a constant and is ignored, and since our final value is relative we don't need to know the unit area)
                  Dr_large <- ((1/nrow(calc.trees[calc.trees$DIA >= 18,]))*sum(calc.trees[calc.trees$DIA >= 18,"DIA"]^Zeide.b))^(1/Zeide.b) #calculate Reineke's diameter for trees greater than 18 inch DBH in condition using Equation 3 of Zeide 1983, aseqation 5 requires an unknown coefficient of variation but achieves the same outcome
                  N_large <- Dr_large^-Zeide.b #use modified equation 1 from Zeide 1983 to arrive at relative density (a is a constant and is ignored, and since our final value is relative we don't need to know the unit area)
                  if(N_large %in% c(Inf,NA,NaN)){N_large <- 0} #if no large trees, we ask the formula to do impossible math. Overcome this limitation by replacing erroneous answer with zero.
                  relative.SDI <- 100*(N_large/N_all) #calculate the proportion of stand density of large trees compared to all trees
                  
                  #assess old growth status based on ecological response unit (ERU)-specific thresholds from Table 9 of MOG document
                  if(ERU %in% c("spruce-fir forest","mixed conifer with aspen","bristlecone pine","pinyon juniper evergreen shrub","pinyon juniper (persistent)","pinyon juniper sagebrush","pinyon juniper deciduous shrub","gambel oak shrubland","arizona walnut","rio grande cottonwood/shrub","narrowleaf cottonwood-spruce, narrowleaf cottonwood/shrub","upper montane conifer/willow")){#begin MOG assessment for (southwest) ERU's that use QMD as a threshold instead of stand density
                    
                    #calculate quadratic mean diameter (QMD) of trees greater than 10 in DBH (adapted from Table 3 of MOG document)
                    #create mature dataframe for several indices
                    sw.mat.df <- local.trees[which(local.trees$DIA >= 10), ] #retain only trees with DBH greater than 10 inch
                    sw.mat.df <- sw.mat.df[which(sw.mat.df$STATUSCD == 1), ] #retain only living trees
                    
                    sw.tpa <- nrow(sw.mat.df)/condition.area #calculate live trees per acre
                    sw.ba <- sum(sw.mat.df$TPA_UNADJ*pi*(sw.mat.df$DIA / 24)*2) #calculate total basal area of trees with DBH greater than 10 in
                    sw.QMD <- sqrt(sw.ba / (sw.tpa * 0.005454)) #calculate quadratic mean diameter of trees with DH greater than 10 in
                    if(sw.QMD %in% c(NA,NaN,Inf)){sw.QMD <- 0} #correct for erroneous answers from impossible math (from lack of trees in sample)
                    
                    local.MOG.status <- ifelse(sw.QMD >= 18, 1, 0) #assign old growth (1) if threshold is met and non-old growth (0) if threshold is not met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end MOG assessment for (southwest) ERU's that use QMD as a threshold instead of stand density
                  
                  if(ERU %in% "mixed conifer - frequent fire"){#begin old growth assessment for (southwest) mixed conifers with frequent fire
                    
                    local.MOG.status <- ifelse(relative.SDI >= 56, 1, 0) #assign old growth (1) if threshold is met and non-old growth (0) if threshold is not met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (southwest) mixed conifers with frequent fire
                  
                  if(ERU %in% c("ponderosa pine forest","ponderosa pine/willow")){#begin old growth assessment for (southwest) ponderosa pine and willow
                    
                    local.MOG.status <- ifelse(relative.SDI >= 57, 1, 0) #assign old growth (1) if threshold is met and non-old growth (0) if threshold is not met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (southwest) ponderosa pine and willow
                  
                  if(ERU %in% c("ponderosa pine - evergreen oak")){#begin old growth assessment for (southwest) ponderosa pine and evergreen oak
                    
                    local.MOG.status <- ifelse(relative.SDI >= 26, 1, 0) #assign old growth (1) if threshold is met and non-old growth (0) if threshold is not met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (southwest) ponderosa pine and evergreen oak
                  
                  if(ERU %in% "pinyon juniper grass"){#begin old growth assessment for (southwest) pinyon juniper grass
                    
                    local.MOG.status <- ifelse(relative.SDI >= 29, 1, 0) #assign old growth (1) if threshold is met and non-old growth (0) if threshold is not met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (southwest) pinyon juniper grass
                  
                  if(ERU %in% c("juniper grass","semi-desert grassland")){#begin old growth assessment for (southwest) juniper grass and semi-desert grassland
                    
                    local.MOG.status <- ifelse(relative.SDI >= 36, 1, 0) #assign old growth (1) if threshold is met and non-old growth (0) if threshold is not met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (southwest) juniper grass and semi-desert grassland
                  
                  if(ERU %in% c("madrean pinyon-oak","madrean encinal woodland")){#begin old growth assessment for (southwest) madrean systems
                    
                    local.MOG.status <- ifelse(relative.SDI >= 20, 1, 0) #assign old growth (1) if threshold is met and non-old growth (0) if threshold is not met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (southwest) madrean systems
                  
                  if(ERU %in% "other"){local.MOG.status <- 0} #if the plot does not fit into an ERU, assign OG status of zero (there is an opportunity to obtain mature status below)
                  
                  #if the plot is not old growth, test to see if it is at least mature based on Table 19 of MOG document
                  if(local.MOG.status == 0){#begin (southwest) mature assessment
                    
                    if(ERU %in% c("arizona walnut","rio grande cottonwood/shrub","gambel oak shrubland","narrowleaf cottonwood-spruce, narrowleaf cottonwood/shrub","upper montane conifer/willow","other") | forest.type %in% c(970:976)){#begin mature assessment for (southwest) hardwoods
                      
                      t.QMDdom <- ifelse(QMDdom >= 3.5, 0.34, 0)
                      t.ddiscore <- ifelse(ddiscore >= 7.7, 0.34, 0)
                      t.HTquart <- ifelse(HTquart >= 10.8, 0.2, 0)
                      t.tpadom <- ifelse(tpadom <= 69.5, 0.12, 0)
                      
                      local.MOG.status <- sum(t.QMDdom,t.ddiscore,t.HTquart,t.tpadom,na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for (southwest) hardwoods
                    
                    if(ERU %in% "juniper grass"){#begin mature assessment for (southwest) juniper grass
                      
                      t.QMDdom <- ifelse(QMDdom >= 10.7, 0.3, 0)
                      t.HTquart <- ifelse(HTquart >= 11.2, 0.27, 0)
                      t.ddiscore <- ifelse(ddiscore >= 19, 0.27, 0)
                      t.HTsd <- ifelse(HTsd >= 4, 0.17, 0)
                      
                      local.MOG.status <- sum(t.QMDdom,t.HTquart,t.ddiscore,t.HTsd,na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for (southwest) juniper grass
                    
                    if(ERU %in% "madrean encinal woodland"){#begin mature assessment for (southwest) madrean encinal woodland
                      
                      t.QMDdom <- ifelse(QMDdom >= 8.8, 0.36, 0)
                      t.HTquart <- ifelse(HTquart >= 15.2, 0.3, 0)
                      t.ddiscore <- ifelse(ddiscore >= 16.8, 0.18, 0)
                      t.tpadom <- ifelse(tpadom <= 56.4, 0.16, 0)
                      
                      local.MOG.status <- sum(t.QMDdom,t.HTquart,t.ddiscore,t.tpadom,na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for (southwest) madrean encinal woodland
                    
                    if(ERU %in% "madrean pinyon-oak"){#begin mature assessment for (southwest) madrean pinyon-oak
                      
                      t.QMDdom <- ifelse(QMDdom >= 8.3, 0.32, 0)
                      t.HTquart <- ifelse(HTquart >= 14.4, 0.28, 0)
                      t.ddiscore <- ifelse(ddiscore >= 23.8, 0.23, 0)
                      t.HTsd <- ifelse(HTsd >= 10.4, 0.16, 0)
                      
                      local.MOG.status <- sum(t.QMDdom,t.HTquart,t.ddiscore,t.HTsd,na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for (southwest) madrean pinyon-oak
                    
                    if(ERU %in% "mixed conifer - frequent fire"){#begin mature assessment for (southwest) mixed conifer - frequent fire
                      
                      t.ddiscore <- ifelse(ddiscore >= 21.4, 0.41, 0)
                      t.QMDdom <- ifelse(QMDdom >= 13.3, 0.28, 0)
                      t.HTsd <- ifelse(HTsd >= 44.7, 0.21, 0)
                      
                      local.MOG.status <- sum(t.ddiscore,t.QMDdom,t.HTsd,na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for (southwest) mixed conifer - frequent fire
                    
                    if(ERU %in% c("mixed conifer with aspen","bristlecone pine")){#begin mature assessment for (southwest) mixed conifer and bristlecone pine
                      
                      t.ddiscore <- ifelse(ddiscore >= 34.5, 0.39, 0)
                      t.HTsd <- ifelse(HTsd >= 41.2, 0.24, 0)
                      t.HTquart <- ifelse(HTquart >= 36.3, 0.22, 0)
                      t.snagbatot <- ifelse(snagbatot <= 15, 0.15, 0)
                      
                      local.MOG.status <- sum(t.ddiscore,t.HTsd,t.HTquart,t.snagbatot,na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for (southwest) mixed conifer and bristlecone pine
                    
                    if(ERU %in% c("pinyon juniper grass","pinyon juniper sagebrush","semi-desert grassland")){#begin mature assessment for (southwest) pinyon juniper grass-sagebursh
                      
                      t.ddiscore <- ifelse(ddiscore >= 19.6, 0.29, 0)
                      t.QMDdom <- ifelse(QMDdom >= 9.5, 0.26, 0)
                      t.HTquart <- ifelse(HTquart >= 12.8, 0.26, 0)
                      t.HTsd <- ifelse(HTsd >= 6.4, 0.19, 0)
                      
                      local.MOG.status <- sum(t.ddiscore,t.QMDdom,t.HTsd,t.HTquart,na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for (southwest)pinyon juniper grass-sagebursh
                    
                    if(ERU %in% c("pinyon juniper (persistent)","pinyon juniper deciduous shrub","pinyon juniper evergreen shrub") | forest.type %in% c(180:185)){#begin mature assessment for (southwest) pinyon shrub - woodland
                      
                      t.ddiscore <- ifelse(ddiscore >= 20.2, 0.46, 0)
                      t.QMDdom <- ifelse(QMDdom >= 9.2, 0.34, 0)
                      t.HTquart <- ifelse(HTquart >= 13.3, 0.21, 0)
                      
                      local.MOG.status <- sum(t.ddiscore,t.QMDdom,t.HTquart,na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for (southwest) pinyon shrub - woodland
                    
                    if(ERU %in% "ponderosa pine forest"){#begin mature assessment for (southwest) ponderosa pine forest
                      
                      t.ddiscore <- ifelse(ddiscore >= 24.3, 0.45, 0)
                      t.badom <- ifelse(badom >= 40, 0.28, 0)
                      t.QMDdom <- ifelse(QMDdom >= 13.5, 0.27, 0)
                      
                      local.MOG.status <- sum(t.ddiscore,t.badom,t.QMDdom,na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for (southwest) ponderosa pine forest
                    
                    if(ERU %in% c("ponderosa pine - evergreen oak","ponderosa pine/willow")){#begin mature assessment for (southwest) ponderosa pine - mixed
                      
                      t.ddiscore <- ifelse(ddiscore >= 32.4, 0.5, 0)
                      t.QMDdom <- ifelse(QMDdom >= 9, 0.32, 0)
                      t.HTsd <- ifelse(HTsd >= 24.1, 0.18, 0)
                      
                      local.MOG.status <- sum(t.ddiscore,t.QMDdom,t.HTsd,na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for (southwest) ponderosa pine - mixed
                    
                    if(ERU %in% "spruce-fir forest" | forest.type %in% c(200:203)){#begin mature assessment for (southwest) spruce/fir
                      
                      t.ddiscore <- ifelse(ddiscore >= 32.4, 0.24, 0)
                      t.HTsd <- ifelse(HTsd >= 51.8, 0.22, 0)
                      t.QMDdom <- ifelse(QMDdom >= 11.4, 0.19, 0)
                      t.HTquart <- ifelse(HTquart >= 43.5, 0.19, 0)
                      t.badom <- ifelse(badom >= 57.4, 0.17, 0)
                      
                      local.MOG.status <- sum(t.ddiscore,t.HTsd,t.QMDdom,t.HTquart,t.badom,na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for (southwest) spruce/fir
                  }#end (southwest) mature assessment
                  
                }#end instructions for southwest region
                
                if(region %in% "pacific southwest"){ #begin instructions for the pacific southwest region
                  
                  #assign regional vegetation type according to Table 12 of MOG document
                  veg.type <- NA #assign as NA to catch errors downstream
                  if(forest.type %in% 341){veg.type <- "coast redwood"}
                  if(forest.type %in% c(371,226,361)){veg.type <- "conifer mixed forests"}
                  if(forest.type %in% 261){veg.type <- "white fir"}
                  if(forest.type %in% c(201,202)){veg.type <- "pacific douglas-fir"}
                  if(forest.type %in% 941){veg.type <- "douglas-fir/tanoak/madrone"}
                  if(forest.type %in% c(241,342,365,366,367)){veg.type <- "mixed subalpine (western white pine assc)"}
                  if(forest.type %in% 270){veg.type <- "mixed subalpine (mountain hemlock assc)"}
                  if(forest.type %in% 369){veg.type <- "mixed subalpine (western juniper assc)"}
                  if(forest.type %in% 901){veg.type <- "mixed subalpine (quaking aspen assc)"}
                  if(forest.type %in% 262){veg.type <- "red fir"}
                  if(forest.type %in% 225){veg.type <- "jeffrey pine"}
                  if(forest.type %in% 281){veg.type <- "lodgepole pine"}
                  if(forest.type %in% 221){veg.type <- "ponderosa pine"}
                  
                  local.trees$tree.age <- ifelse(!(is.na(local.trees$BHAGE)),local.trees$BHAGE, #assess tree age (could be stored in several locations)
                                                 ifelse(!(is.na(local.trees$TOTAGE)),local.trees$TOTAGE,raw.stand.age)) 
                  site.index <- local.trees$HT * (0.25489 + (29.377 / local.trees$tree.age))
                  site.index <- mean(site.index, na.rm=TRUE) #calculate the mean site index
                  if(!is.na(site.index)){site.index <- ifelse(site.index < 45, "low","productive")} #assign to binary growth class based on Dunning's site index cutoff of 45
                  
                  
                  if(is.na(site.index)){#begin secondary site index classification if no values were available to calculate site index
                    
                    #need to get site tree to condition level and look at crown class
                    site.clcd <- max(local.condition$SITECLCD,local.condition$SITECLCDEST,na.rm=TRUE) #use either the measured or estimated value of productivity (one of these should always be NA, or else they should be the same value)
                    
                    if(!(is.na(site.clcd))){site.index <- ifelse(site.clcd >= 5, "low", "productive")}
                    
                  }#end secondary site index classification if no values were available to calculate site index
                  
                  
                  #determine old growth status per Table 12 of the MOG document
                  if(veg.type %in% "coast redwood"){#begin old growth assessment for (pacific southwest) coastal redwood
                    
                    #there is no productivity requirement for this forest type
                    
                    local.large.trees <- local.trees[local.trees$DIA >= 40, ] #isolate to large trees using threshold specific to the forest type
                    large.trees.p.acre <- nrow(local.large.trees)/condition.area #calculate large trees per acre
                    
                    #there is no age requirement for this forest type
                    
                    local.MOG.status <- ifelse(large.trees.p.acre >= 15, 1, 0) #assign old growth status (1) if threshold is met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (pacific southwest) coastal redwood
                  
                  if(veg.type %in% "conifer mixed forests" & site.index %in% "productive"){#begin old growth assessment for (pacific southwest) productive conifer mixed forests
                    
                    local.large.trees <- local.trees[local.trees$DIA >= 39, ] #isolate to large trees using threshold specific to the forest type
                    large.trees.p.acre <- nrow(local.large.trees)/condition.area #calculate large trees per acre
                    large.trees.p.acre <- ifelse(large.trees.p.acre >= 6, 1, 0)
                    
                    max.age <- max(c(raw.stand.age,local.trees$tree.age), na.rm=TRUE) #determine the age of the oldest tree in the stand
                    max.age <- ifelse(max.age >= 188, 1, 0)
                    
                    local.MOG.status <- ifelse(sum(large.trees.p.acre, max.age, na.rm=TRUE) == 2, 1, 0) #assign old growth status (1) if threshold is met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (pacific southwest) productive conifer mixed forests
                  
                  if(veg.type %in% "conifer mixed forests" & site.index %in% "low"){#begin old growth assessment for (pacific southwest) low conifer mixed forests
                    
                    local.large.trees <- local.trees[local.trees$DIA >= 29, ] #isolate to large trees using threshold specific to the forest type
                    large.trees.p.acre <- nrow(local.large.trees)/condition.area #calculate large trees per acre
                    large.trees.p.acre <- ifelse(large.trees.p.acre >= 5, 1, 0)
                    
                    max.age <- max(c(raw.stand.age,local.trees$tree.age), na.rm=TRUE) #determine the age of the oldest tree in the stand
                    max.age <- ifelse(max.age >= 256, 1, 0)
                    
                    local.MOG.status <- ifelse(sum(large.trees.p.acre, max.age, na.rm=TRUE) == 2, 1, 0) #assign old growth status (1) if threshold is met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (pacific southwest) low conifer mixed forests
                  
                  if(veg.type %in% "white fir"){#begin old growth assessment for (pacific southwest) white fir
                    
                    #determine if the plot falls within the northwest forest plan (this applies to white fir only)
                    nwfp.test <- st_filter(x = local.plot,
                                           y = nwfp)
                    nwfp.test <- ifelse(nrow(nwfp.test) > 0, "nwfp", "non-nwfp")
                    
                    if(nwfp.test %in% "nwfp"){#begin sub-assessment for (pacific southwest) white fir within the northwest forest plan boundary
                      
                      if(site.index %in% "productive"){#begin sub-sub-assessment of (pacific southwest) productive white fir within the nwfp boundary
                        
                        local.large.trees <- local.trees[local.trees$DIA >= 30, ] #isolate to large trees using threshold specific to the forest type
                        large.trees.p.acre <- nrow(local.large.trees)/condition.area #calculate large trees per acre
                        large.trees.p.acre <- ifelse(large.trees.p.acre >= 5, 1, 0)
                        
                        max.age <- max(c(raw.stand.age,local.trees$tree.age), na.rm=TRUE) #determine the age of the oldest tree in the stand
                        max.age <- ifelse(max.age >= 160, 1, 0)
                        
                        local.MOG.status <- ifelse(sum(large.trees.p.acre, max.age, na.rm=TRUE) == 2, 1, 0) #assign old growth status (1) if threshold is met
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#begin sub-sub-assessment of (pacific southwest) productive white fir within the nwfp boundary
                      
                      if(site.index %in% "low"){#begin sub-sub-assessment of (pacific southwest) low white fir within the nwfp boundary
                        
                        local.large.trees <- local.trees[local.trees$DIA >= 25, ] #isolate to large trees using threshold specific to the forest type
                        large.trees.p.acre <- nrow(local.large.trees)/condition.area #calculate large trees per acre
                        large.trees.p.acre <- ifelse(large.trees.p.acre >= 23, 1, 0)
                        
                        max.age <- max(c(raw.stand.age,local.trees$tree.age), na.rm=TRUE) #determine the age of the oldest tree in the stand
                        max.age <- ifelse(max.age >= 303, 1, 0)
                        
                        local.MOG.status <- ifelse(sum(large.trees.p.acre, max.age, na.rm=TRUE) == 2, 1, 0) #assign old growth status (1) if threshold is met
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#begin sub-sub-assessment of (pacific southwest) low white fir within the nwfp boundary
                      
                    }#end sub-assessment for (pacific southwest) white fir within the northwest forest plan boundary
                    
                    if(nwfp.test %in% "non-nwfp"){#begin sub-assessment for (pacific southwest) white fir outside the northwest forest plan boundary
                      
                      if(site.index %in% "productive"){#begin sub-sub-assessment of (pacific southwest) productive white fir outside the nwfp boundary
                        
                        local.large.trees <- local.trees[local.trees$DIA >= 39, ] #isolate to large trees using threshold specific to the forest type
                        large.trees.p.acre <- nrow(local.large.trees)/condition.area #calculate large trees per acre
                        large.trees.p.acre <- ifelse(large.trees.p.acre >= 6, 1, 0)
                        
                        max.age <- max(c(raw.stand.age,local.trees$tree.age), na.rm=TRUE) #determine the age of the oldest tree in the stand
                        max.age <- ifelse(max.age >= 143, 1, 0)
                        
                        local.MOG.status <- ifelse(sum(large.trees.p.acre, max.age, na.rm=TRUE) == 2, 1, 0) #assign old growth status (1) if threshold is met
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end sub-sub-assessment of (pacific southwest) productive white fir outside the nwfp boundary
                      
                      if(site.index %in% "low"){#begin sub-sub-assessment of (pacific southwest) low white fir outside the nwfp boundary
                        
                        local.large.trees <- local.trees[local.trees$DIA >= 29, ] #isolate to large trees using threshold specific to the forest type
                        large.trees.p.acre <- nrow(local.large.trees)/condition.area #calculate large trees per acre
                        large.trees.p.acre <- ifelse(large.trees.p.acre >= 8, 1, 0)
                        
                        max.age <- max(c(raw.stand.age,local.trees$tree.age), na.rm=TRUE) #determine the age of the oldest tree in the stand
                        max.age <- ifelse(max.age >= 239, 1, 0)
                        
                        local.MOG.status <- ifelse(sum(large.trees.p.acre, max.age, na.rm=TRUE) == 2, 1, 0) #assign old growth status (1) if threshold is met
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end sub-sub-assessment of (pacific southwest) low white fir outside the nwfp boundary
                    }#end sub-assessment for (pacific southwest) white fir within the northwest forest plan boundary
                  }#end old growth assessment for (pacific southwest) white fir
                  
                  if(veg.type %in% "pacific douglas-fir" & site.index %in% "productive"){#begin old growth assessment for (pacific southwest) productive pacific douglas-fir
                    
                    local.large.trees <- local.trees[local.trees$DIA >= 40, ] #isolate to large trees using threshold specific to the forest type
                    large.trees.p.acre <- nrow(local.large.trees)/condition.area #calculate large trees per acre
                    large.trees.p.acre <- ifelse(large.trees.p.acre >= 12, 1, 0)
                    
                    max.age <- max(c(raw.stand.age,local.trees$tree.age), na.rm=TRUE) #determine the age of the oldest tree in the stand
                    max.age <- ifelse(max.age >= 180, 1, 0)
                    
                    local.MOG.status <- ifelse(sum(large.trees.p.acre, max.age, na.rm=TRUE) == 2, 1, 0) #assign old growth status (1) if threshold is met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (pacific southwest) productive pacific douglas-fir
                  
                  if(veg.type %in% "pacific douglas-fir" & site.index %in% "low"){#begin old growth assessment for (pacific southwest) low pacific douglas-fir
                    
                    local.large.trees <- local.trees[local.trees$DIA >= 30, ] #isolate to large trees using threshold specific to the forest type
                    large.trees.p.acre <- nrow(local.large.trees)/condition.area #calculate large trees per acre
                    large.trees.p.acre <- ifelse(large.trees.p.acre >= 18, 1, 0)
                    
                    max.age <- max(c(raw.stand.age,local.trees$tree.age), na.rm=TRUE) #determine the age of the oldest tree in the stand
                    max.age <- ifelse(max.age >= 260, 1, 0)
                    
                    local.MOG.status <- ifelse(sum(large.trees.p.acre, max.age, na.rm=TRUE) == 2, 1, 0) #assign old growth status (1) if threshold is met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (pacific southwest) low pacific douglas-fir
                  
                  if(veg.type %in% "douglas-fir/tanoak/madrone" & site.index %in% "productive"){#begin old growth assessment for (pacific southwest) productive douglas-fir/tanoak/madrone
                    
                    local.large.trees <- local.trees[local.trees$DIA >= 30, ] #isolate to large trees using threshold specific to the forest type
                    large.trees.p.acre <- nrow(local.large.trees)/condition.area #calculate large trees per acre
                    large.trees.p.acre <- ifelse(large.trees.p.acre >= 10, 1, 0)
                    
                    max.age <- max(c(raw.stand.age,local.trees$tree.age), na.rm=TRUE) #determine the age of the oldest tree in the stand
                    max.age <- ifelse(max.age >= 180, 1, 0)
                    
                    local.MOG.status <- ifelse(sum(large.trees.p.acre, max.age, na.rm=TRUE) == 2, 1, 0) #assign old growth status (1) if threshold is met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (pacific southwest) productive douglas-fir/tanoak/madrone
                  
                  if(veg.type %in% "douglas-fir/tanoak/madrone" & site.index %in% "low"){#begin old growth assessment for (pacific southwest) low douglas-fir/tanoak/madrone
                    
                    local.large.trees <- local.trees[local.trees$DIA >= 30, ] #isolate to large trees using threshold specific to the forest type
                    large.trees.p.acre <- nrow(local.large.trees)/condition.area #calculate large trees per acre
                    large.trees.p.acre <- ifelse(large.trees.p.acre >= 8, 1, 0)
                    
                    max.age <- max(c(raw.stand.age,local.trees$tree.age), na.rm=TRUE) #determine the age of the oldest tree in the stand
                    max.age <- ifelse(max.age >= 300, 1, 0)
                    
                    local.MOG.status <- ifelse(sum(large.trees.p.acre, max.age, na.rm=TRUE) == 2, 1, 0) #assign old growth status (1) if threshold is met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (pacific southwest) low douglas-fir/tanoak/madrone
                  
                  if(veg.type %in% "mixed subalpine (western white pine assc)" & site.index %in% "productive"){#begin old growth assessment for (pacific southwest) productive mixed subalpine (western white pine assc)
                    
                    local.large.trees <- local.trees[local.trees$DIA >= 30, ] #isolate to large trees using threshold specific to the forest type
                    large.trees.p.acre <- nrow(local.large.trees)/condition.area #calculate large trees per acre
                    large.trees.p.acre <- ifelse(large.trees.p.acre >= 9, 1, 0)
                    
                    max.age <- max(c(raw.stand.age,local.trees$tree.age), na.rm=TRUE) #determine the age of the oldest tree in the stand
                    max.age <- ifelse(max.age >= 150, 1, 0)
                    
                    local.MOG.status <- ifelse(sum(large.trees.p.acre, max.age, na.rm=TRUE) == 2, 1, 0) #assign old growth status (1) if threshold is met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (pacific southwest) productive mixed subalpine (western white pine assc)
                  
                  if(veg.type %in% "mixed subalpine (western white pine assc)" & site.index %in% "low"){#begin old growth assessment for (pacific southwest) low mixed subalpine (western white pine assc)
                    
                    local.large.trees <- local.trees[local.trees$DIA >= 30, ] #isolate to large trees using threshold specific to the forest type
                    large.trees.p.acre <- nrow(local.large.trees)/condition.area #calculate large trees per acre
                    large.trees.p.acre <- ifelse(large.trees.p.acre >= 10, 1, 0)
                    
                    max.age <- max(c(raw.stand.age,local.trees$tree.age), na.rm=TRUE) #determine the age of the oldest tree in the stand
                    max.age <- ifelse(max.age >= 200, 1, 0)
                    
                    local.MOG.status <- ifelse(sum(large.trees.p.acre, max.age, na.rm=TRUE) == 2, 1, 0) #assign old growth status (1) if threshold is met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (pacific southwest) low mixed subalpine (western white pine assc)
                  
                  if(veg.type %in% "mixed subalpine (mountain hemlock assc)" & site.index %in% "productive"){#begin old growth assessment for (pacific southwest) productive mixed subalpine (mountain hemlock assc)
                    
                    local.large.trees <- local.trees[local.trees$DIA >= 30, ] #isolate to large trees using threshold specific to the forest type
                    large.trees.p.acre <- nrow(local.large.trees)/condition.area #calculate large trees per acre
                    large.trees.p.acre <- ifelse(large.trees.p.acre >= 12, 1, 0)
                    
                    max.age <- max(c(raw.stand.age,local.trees$tree.age), na.rm=TRUE) #determine the age of the oldest tree in the stand
                    max.age <- ifelse(max.age >= 150, 1, 0)
                    
                    local.MOG.status <- ifelse(sum(large.trees.p.acre, max.age, na.rm=TRUE) == 2, 1, 0) #assign old growth status (1) if threshold is met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (pacific southwest) productive mixed subalpine (mountain hemlock assc)
                  
                  if(veg.type %in% "mixed subalpine (mountain hemlock assc)" & site.index %in% "low"){#begin old growth assessment for (pacific southwest) low mixed subalpine (mountain hemlock assc)
                    
                    local.large.trees <- local.trees[local.trees$DIA >= 30, ] #isolate to large trees using threshold specific to the forest type
                    large.trees.p.acre <- nrow(local.large.trees)/condition.area #calculate large trees per acre
                    large.trees.p.acre <- ifelse(large.trees.p.acre >= 6, 1, 0)
                    
                    max.age <- max(c(raw.stand.age,local.trees$tree.age), na.rm=TRUE) #determine the age of the oldest tree in the stand
                    max.age <- ifelse(max.age >= 200, 1, 0)
                    
                    local.MOG.status <- ifelse(sum(large.trees.p.acre, max.age, na.rm=TRUE) == 2, 1, 0) #assign old growth status (1) if threshold is met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (pacific southwest) low mixed subalpine (mountain hemlock assc)
                  
                  if(veg.type %in% "mixed subalpine (western juniper assc)"){#begin old growth assessment for (pacific southwest) mixed subalpine (western juniper assc)
                    
                    #there is no productivity requirement for this forest type
                    
                    local.large.trees <- local.trees[local.trees$DIA >= 30, ] #isolate to large trees using threshold specific to the forest type
                    large.trees.p.acre <- nrow(local.large.trees)/condition.area #calculate large trees per acre
                    large.trees.p.acre <- ifelse(large.trees.p.acre >= 5, 1, 0)
                    
                    max.age <- max(c(raw.stand.age,local.trees$tree.age), na.rm=TRUE) #determine the age of the oldest tree in the stand
                    max.age <- ifelse(max.age >= 200, 1, 0)
                    
                    local.MOG.status <- ifelse(sum(large.trees.p.acre, max.age, na.rm=TRUE) == 2, 1, 0) #assign old growth status (1) if threshold is met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (pacific southwest) mixed subalpine (western juniper assc)
                  
                  if(veg.type %in% "mixed subalpine (quaking aspen assc)" & site.index %in% "productive"){#begin old growth assessment for (pacific southwest) productive mixed subalpine (quaking aspen assc)
                    
                    local.large.aspens <- local.trees[local.trees$SPGRPCD %in% 44 & local.trees$DIA >= 18, ] #isolate to large trees using threshold specific to the forest type, taking care to use different diameter thresholds for different species
                    local.large.conifer <- local.trees[local.trees$SPGRPCD %in% c(1:24) & local.trees$DIA >= 30, ] #isolate to large trees using threshold specific to the forest type, taking care to use different diameter thresholds for different species
                    local.large.trees <- rbind(local.large.aspens,local.large.conifer)
                    large.trees.p.acre <- nrow(local.large.trees)/condition.area #calculate large trees per acre
                    large.trees.p.acre <- ifelse(large.trees.p.acre >= 5, 1, 0)
                    
                    max.age <- max(c(raw.stand.age,local.trees$tree.age), na.rm=TRUE) #determine the age of the oldest tree in the stand
                    max.age <- ifelse(max.age >= 80, 1, 0)
                    
                    local.MOG.status <- ifelse(sum(large.trees.p.acre, max.age, na.rm=TRUE) == 2, 1, 0) #assign old growth status (1) if threshold is met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (pacific southwest) productive mixed subalpine (quaking aspen assc)
                  
                  if(veg.type %in% "mixed subalpine (quaking aspen assc)" & site.index %in% "low"){#begin old growth assessment for (pacific southwest) low mixed subalpine (quaking aspen assc)
                    
                    local.large.aspens <- local.trees[local.trees$SPGRPCD %in% 44 & local.trees$DIA >= 18, ] #isolate to large trees using threshold specific to the forest type, taking care to use different diameter thresholds for different species
                    local.large.conifer <- local.trees[local.trees$SPGRPCD %in% c(1:24) & local.trees$DIA >= 30, ] #isolate to large trees using threshold specific to the forest type, taking care to use different diameter thresholds for different species
                    local.large.trees <- rbind(local.large.aspens,local.large.conifer)
                    large.trees.p.acre <- nrow(local.large.trees)/condition.area #calculate large trees per acre
                    large.trees.p.acre <- ifelse(large.trees.p.acre >= 1, 1, 0)
                    
                    max.age <- max(c(raw.stand.age,local.trees$tree.age), na.rm=TRUE) #determine the age of the oldest tree in the stand
                    max.age <- ifelse(max.age >= 80, 1, 0)
                    
                    local.MOG.status <- ifelse(sum(large.trees.p.acre, max.age, na.rm=TRUE) == 2, 1, 0) #assign old growth status (1) if threshold is met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (pacific southwest) low mixed subalpine (quaking aspen assc)
                  
                  if(veg.type %in% "red fir" & site.index %in% "productive"){#begin old growth assessment for (pacific southwest) productive red fir
                    
                    local.large.trees <- local.trees[local.trees$DIA >= 30, ] #isolate to large trees using threshold specific to the forest type
                    large.trees.p.acre <- nrow(local.large.trees)/condition.area #calculate large trees per acre
                    large.trees.p.acre <- ifelse(large.trees.p.acre >= 8, 1, 0)
                    
                    max.age <- max(c(raw.stand.age,local.trees$tree.age), na.rm=TRUE) #determine the age of the oldest tree in the stand
                    max.age <- ifelse(max.age >= 150, 1, 0)
                    
                    local.MOG.status <- ifelse(sum(large.trees.p.acre, max.age, na.rm=TRUE) == 2, 1, 0) #assign old growth status (1) if threshold is met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (pacific southwest) productive red fir
                  
                  if(veg.type %in% "red fir" & site.index %in% "low"){#begin old growth assessment for (pacific southwest) low red fir
                    
                    local.large.trees <- local.trees[local.trees$DIA >= 36, ] #isolate to large trees using threshold specific to the forest type
                    large.trees.p.acre <- nrow(local.large.trees)/condition.area #calculate large trees per acre
                    large.trees.p.acre <- ifelse(large.trees.p.acre >= 5, 1, 0)
                    
                    max.age <- max(c(raw.stand.age,local.trees$tree.age), na.rm=TRUE) #determine the age of the oldest tree in the stand
                    max.age <- ifelse(max.age >= 200, 1, 0)
                    
                    local.MOG.status <- ifelse(sum(large.trees.p.acre, max.age, na.rm=TRUE) == 2, 1, 0) #assign old growth status (1) if threshold is met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (pacific southwest) low red fir
                  
                  if(veg.type %in% "jeffrey pine" & site.index %in% "productive"){#begin old growth assessment for (pacific southwest) productive jeffrey pine
                    
                    local.large.trees <- local.trees[local.trees$DIA >= 30, ] #isolate to large trees using threshold specific to the forest type
                    large.trees.p.acre <- nrow(local.large.trees)/condition.area #calculate large trees per acre
                    large.trees.p.acre <- ifelse(large.trees.p.acre >= 3, 1, 0)
                    
                    max.age <- max(c(raw.stand.age,local.trees$tree.age), na.rm=TRUE) #determine the age of the oldest tree in the stand
                    max.age <- ifelse(max.age >= 150, 1, 0)
                    
                    local.MOG.status <- ifelse(sum(large.trees.p.acre, max.age, na.rm=TRUE) == 2, 1, 0) #assign old growth status (1) if threshold is met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (pacific southwest) productive jeffrey pine
                  
                  if(veg.type %in% "jeffrey pine" & site.index %in% "low"){#begin old growth assessment for (pacific southwest) low jeffrey pine
                    
                    local.large.trees <- local.trees[local.trees$DIA >= 30, ] #isolate to large trees using threshold specific to the forest type
                    large.trees.p.acre <- nrow(local.large.trees)/condition.area #calculate large trees per acre
                    large.trees.p.acre <- ifelse(large.trees.p.acre >= 1, 1, 0)
                    
                    max.age <- max(c(raw.stand.age,local.trees$tree.age), na.rm=TRUE) #determine the age of the oldest tree in the stand
                    max.age <- ifelse(max.age >= 200, 1, 0)
                    
                    local.MOG.status <- ifelse(sum(large.trees.p.acre, max.age, na.rm=TRUE) == 2, 1, 0) #assign old growth status (1) if threshold is met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (pacific southwest) low jeffrey pine
                  
                  if(veg.type %in% "lodgepole pine" & site.index %in% "productive"){#begin old growth assessment for (pacific southwest) productive lodgepole pine
                    
                    local.large.trees <- local.trees[local.trees$DIA >= 36, ] #isolate to large trees using threshold specific to the forest type
                    large.trees.p.acre <- nrow(local.large.trees)/condition.area #calculate large trees per acre
                    large.trees.p.acre <- ifelse(large.trees.p.acre >= 7, 1, 0)
                    
                    max.age <- max(c(raw.stand.age,local.trees$tree.age), na.rm=TRUE) #determine the age of the oldest tree in the stand
                    max.age <- ifelse(max.age >= 150, 1, 0)
                    
                    local.MOG.status <- ifelse(sum(large.trees.p.acre, max.age, na.rm=TRUE) == 2, 1, 0) #assign old growth status (1) if threshold is met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (pacific southwest) productive lodgepole pine
                  
                  if(veg.type %in% "lodgepole pine" & site.index %in% "low"){#begin old growth assessment for (pacific southwest) low lodgepole pine
                    
                    local.large.trees <- local.trees[local.trees$DIA >= 36, ] #isolate to large trees using threshold specific to the forest type
                    large.trees.p.acre <- nrow(local.large.trees)/condition.area #calculate large trees per acre
                    large.trees.p.acre <- ifelse(large.trees.p.acre >= 4, 1, 0)
                    
                    max.age <- max(c(raw.stand.age,local.trees$tree.age), na.rm=TRUE) #determine the age of the oldest tree in the stand
                    max.age <- ifelse(max.age >= 200, 1, 0)
                    
                    local.MOG.status <- ifelse(sum(large.trees.p.acre, max.age, na.rm=TRUE) == 2, 1, 0) #assign old growth status (1) if threshold is met
                    MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    
                  }#end old growth assessment for (pacific southwest) low lodgepole pine
                  
                  pond.class <- NA #assigning an empty ponderosa class to catch downstream errors
                  if(veg.type %in% "ponderosa pine"){#begin old growth assessment for (pacific southwest) ponderosa pine
                    
                    ECOSUBCD <- FIA.geomplot[which(FIA.geomplot$CN %in% local.plot$CN),"ECOSUBCD"] #assign ecological subregion
                    pond.class <- ifelse(ECOSUBCD %in% c("M261Ea","M261Eb","M261Ec","M261Ei","M261Ej"), "interior", #determine if subregion is interior or pacific classification per footnote A in Table 12 of MOG document
                                         ifelse(substr(ECOSUBCD, 1,5) %in% c("M261G","M261D"), "interior",
                                                ifelse(substr(ECOSUBCD, 1,4) %in% "342B", "interior", "pacific")))
                    if(ECOSUBCD %in% "M261Di,M"){pond.class <- "pacific"} #reassign based on a particular exemption (see Footnote A in Tbale 12 of MOG document)
                    
                    if(pond.class %in% "interior" & site.index %in% "productive"){ #begin subassessment for (pacific southwest) productive interior ponderosa pine
                      local.large.trees <- local.trees[local.trees$DIA >= 21, ] #isolate to large trees using threshold specific to the forest type
                      large.trees.p.acre <- nrow(local.large.trees)/condition.area #calculate large trees per acre
                      large.trees.p.acre <- ifelse(large.trees.p.acre >= 19, 1, 0)
                      
                      max.age <- max(c(raw.stand.age,local.trees$tree.age), na.rm=TRUE) #determine the age of the oldest tree in the stand
                      max.age <- ifelse(max.age >= 150, 1, 0)
                      
                      local.MOG.status <- ifelse(sum(large.trees.p.acre, max.age, na.rm=TRUE) == 2, 1, 0) #assign old growth status (1) if threshold is met
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    } #end subassessment for (pacific southwest) productive interior ponderosa pine
                    
                    if(pond.class %in% "interior" & site.index %in% "low"){ #begin subassessment for (pacific southwest) low interior ponderosa pine
                      local.large.trees <- local.trees[local.trees$DIA >= 21, ] #isolate to large trees using threshold specific to the forest type
                      large.trees.p.acre <- nrow(local.large.trees)/condition.area #calculate large trees per acre
                      large.trees.p.acre <- ifelse(large.trees.p.acre >= 16, 1, 0)
                      
                      max.age <- max(c(raw.stand.age,local.trees$tree.age), na.rm=TRUE) #determine the age of the oldest tree in the stand
                      max.age <- ifelse(max.age >= 200, 1, 0)
                      
                      local.MOG.status <- ifelse(sum(large.trees.p.acre, max.age, na.rm=TRUE) == 2, 1, 0) #assign old growth status (1) if threshold is met
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    } #end subassessment for (pacific southwest) low interior ponderosa pine
                    
                    if(pond.class %in% "pacific"){ #begin subassessment for (pacific southwest) pacific ponderosa pine
                      
                      #there is no productivity qualification for pacific ponderosas
                      
                      local.large.trees <- local.trees[local.trees$DIA >= 30, ] #isolate to large trees using threshold specific to the forest type
                      large.trees.p.acre <- nrow(local.large.trees)/condition.area #calculate large trees per acre
                      large.trees.p.acre <- ifelse(large.trees.p.acre >= 9, 1, 0)
                      
                      max.age <- max(c(raw.stand.age,local.trees$tree.age), na.rm=TRUE) #determine the age of the oldest tree in the stand
                      max.age <- ifelse(max.age >= 125, 1, 0)
                      
                      local.MOG.status <- ifelse(sum(large.trees.p.acre, max.age, na.rm=TRUE) == 2, 1, 0) #assign old growth status (1) if threshold is met
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                    } #end subassessment for (pacific southwest) pacific ponderosa pine
                  }#end old growth assessment for (pacific southwest) ponderosa pine
                  
                  #if the plot is not old growth, test to see if it is at least mature based on Table 19 of MOG document
                  if(local.MOG.status == 0){#begin (pacific southwest) mature assessment
                    
                    if(veg.type %in% "douglas-fir/tanoak/madrone"){#begin mature assessment for (pacific southwest) douglas-fir/tanoak/madrone
                      
                      t.ddiscore <- ifelse(ddiscore >= 53.3, 0.45, 0)
                      t.QMDdom <- ifelse(QMDdom >= 14.8, 0.29, 0)
                      t.tpadom <- ifelse(tpadom <= 76.6, 0.25, 0)
                      
                      local.MOG.status <- sum(t.ddiscore,t.QMDdom,t.tpadom,na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for (pacific southwest) douglas-fir/tanoak/madrone
                    
                    if(veg.type %in% "jeffrey pine"){#begin mature assessment for (pacific southwest) jeffery pine
                      
                      t.QMDdom <- ifelse(QMDdom >= 10.3, 0.52, 0)
                      t.ddiscore <- ifelse(ddiscore >= 30.8, 0.25, 0)
                      t.HTsd <- ifelse(HTsd >= 31.5, 0.23, 0)
                      
                      local.MOG.status <- sum(t.ddiscore,t.QMDdom,t.HTsd,na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for (pacific southwest) jeffery pine
                    
                    if(veg.type %in% c("conifer mixed forests","lodgepole pine","mixed subalpine (western white pine assc)","mixed subalpine (mountain hemlock assc)")){#begin mature assessment for (pacific southwest) mixed conifer
                      
                      t.QMDdom <- ifelse(QMDdom >= 13.1, 0.6, 0)
                      t.ddiscore <- ifelse(ddiscore >= 42.1, 0.4, 0)
                      
                      local.MOG.status <- sum(t.ddiscore,t.QMDdom,na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for (pacific southwest) mixed conifer
                    
                    if(veg.type %in% "ponderosa pine" & pond.class %in% "interior"){#begin mature assessment for a special case of (pacific southwest) mixed conifer
                      
                      t.QMDdom <- ifelse(QMDdom >= 13.1, 0.6, 0)
                      t.ddiscore <- ifelse(ddiscore >= 42.1, 0.4, 0)
                      
                      local.MOG.status <- sum(t.ddiscore,t.QMDdom,na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for a special case of (pacific southwest) mixed conifer
                    
                    if(veg.type %in% c("coast redwood","pacific douglas-fir")){#begin mature assessment for (pacific southwest) pacific conifer
                      
                      t.ddiscore <- ifelse(ddiscore >= 52.6, 0.4, 0)
                      t.QMDdom <- ifelse(QMDdom >= 25.3, 0.35, 0)
                      t.snagbatot <- ifelse(snagbatot >= 2.7, 0.26, 0)
                      
                      local.MOG.status <- sum(t.ddiscore,t.QMDdom,t.snagbatot,na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for (pacific southwest) pacific conifer
                    
                    if(veg.type %in% "ponderosa pine" & pond.class %in% "pacific"){#begin mature assessment for a special case of (pacific southwest) pacific conifer
                      
                      t.ddiscore <- ifelse(ddiscore >= 52.6, 0.4, 0)
                      t.QMDdom <- ifelse(QMDdom >= 25.3, 0.35, 0)
                      t.snagbatot <- ifelse(snagbatot >= 2.7, 0.26, 0)
                      
                      local.MOG.status <- sum(t.ddiscore,t.QMDdom,t.snagbatot,na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for a special case of (pacific southwest) pacific conifer
                    
                    if(veg.type %in% "red fir"){#begin mature assessment for (pacific southwest) red fir
                      
                      t.ddiscore <- ifelse(ddiscore >= 48.3, 0.32, 0)
                      t.QMDdom <- ifelse(QMDdom >= 18.1, 0.28, 0)
                      t.HTquart <- ifelse(HTquart >= 66.2, 0.23, 0)
                      t.HTsd <- ifelse(HTsd >= 43.6, 0.17, 0)
                      
                      local.MOG.status <- sum(t.ddiscore,t.QMDdom,t.HTquart,t.HTsd,na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for (pacific southwest) red fir
                    
                    if(veg.type %in% "white fir"){#begin mature assessment for (pacific southwest) white fir
                      
                      t.ddiscore <- ifelse(ddiscore >= 47.5, 0.31, 0)
                      t.HTquart <- ifelse(HTquart >= 68.5, 0.31, 0)
                      t.badom <- ifelse(badom >= 150, 0.21, 0)
                      t.snagbatot <- ifelse(snagbatot >= 24.9, 0.16, 0)
                      
                      local.MOG.status <- sum(t.ddiscore,t.HTquart,t.badom,t.snagbatot,na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for (pacific southwest) white fir
                    
                    if(veg.type %in% c("mixed subalpine (quaking aspen assc)") | forest.type %in% c(910:912, 940:943, 700:722, 920:935, 960:976, 900:905)){#begin mature assessment for (pacific southwest) hardwoods
                      
                      t.ddiscore <- ifelse(ddiscore >= 47.5, 0.31, 0)
                      t.HTquart <- ifelse(HTquart >= 68.5, 0.31, 0)
                      t.badom <- ifelse(badom >= 150, 0.21, 0)
                      t.snagbatot <- ifelse(snagbatot >= 24.9, 0.16, 0)
                      
                      local.MOG.status <- sum(t.ddiscore,t.HTquart,t.badom,t.snagbatot,na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for (pacific southwest) hardwoods
                    
                    if(veg.type %in% c("mixed subalpine (western juniper assc)") | forest.type %in% c(180:185, 360:369)){#begin mature assessment for (pacific southwest) softwoods
                      
                      t.QMDdom <- ifelse(QMDdom >= 14.2, 0.54, 0)
                      t.badom <- ifelse(badom >= 30.9, 0.46, 0)
                      
                      local.MOG.status <- sum(t.QMDdom,t.badom,na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for (pacific southwest) softwoods
                  }#end (pacific southwest) mature assessment
                  
                  #conditions[c,"MOG"] <- max(MOG.vector,na.rm=TRUE) #if MOG conditions are met for any of the forest types assigned to the plot, assign 1 for OG. Otherwise assign 0.5 for mature or zero for neither
                  #conditions[c,"age"] <- raw.stand.age #save stand age
                  
                } #end instructions for the pacific southwest region
                
                if(region %in% "pacific northwest"){ #begin instructions for pacific northwest region
                  
                  #determine plant association zone
                  temp.FIA.plot <- project(x = vect(local.plot), #transform local plot to match crs of plant association zone (this is faster than trying to transform PAZ)
                                           y = paz)
                  paz.val <- terra::extract(y = temp.FIA.plot, #extract value of cell overlapping the point
                                            x = paz)
                  paz.val <- paz.val[1,2] #reduce to just the value
                  
                  paz.group <- NA #assign NA to accommodate downstream errors
                  if(paz.val %in% c(56,57,58,59)){paz.group <- "white fir - grand fir"}
                  if(paz.val %in% c(40,41,42,43)){paz.group <- "douglas fir"}
                  if(paz.val %in% c(25,27,29)){paz.group <- "lodgepole pine"}
                  if(paz.val %in% c(87,88,89)){paz.group <- "silver fir"}
                  if(paz.val %in% c(30,31,33,34)){paz.group <- "ponderosa pine"}
                  if(paz.val %in% c(61,62,63,64)){paz.group <- "subalpine fir"}
                  if(paz.val %in% c(82,83,84)){paz.group <- "western hemlock"}
                  if(paz.val %in% c(91,92,93,94)){paz.group <- "mountain hemlock"}
                  if(paz.val %in% c(21,22)){paz.group <- "juniper"}
                  if(paz.val %in% c(18)){paz.group <- "oak woodland"}
                  if(paz.val %in% c(79)){paz.group <- "port orford cedar"}
                  if(paz.val %in% c(48,49)){paz.group <- "redwood"}
                  if(paz.val %in% c(71,72,73)){paz.group <- "shasta red fir"}
                  if(paz.val %in% c(46,47)){paz.group <- "sitka spruce"}
                  if(paz.val %in% c(51,52)){paz.group <- "tanoak"}
                  if(paz.val %in% c(37,38,39)){paz.group <- "jeffrey pine - knobcone pine"}
                  
                  
                  nwfp.test <- st_filter(x = local.plot, #test if plot resides within the northwest forest plan
                                         y = nwfp)
                  
                  #begin (pacific northwest) old growth assessment for areas INSIDE of the northwest forest. This is documented in Table 13 of MOG document, but instructions are lost in translation. The following process was based on personal communication with regional staff familiar with the process.
                  if(nrow(nwfp.test) > 0){#begin instructions for if the plot is inside the northwest forest plan
                    
                    #test for at least 10% live tree cover (due to forest fires)
                    local.live.trees <- local.trees[local.trees$STATUSCD == 1, ] #subset to live trees only
                    p.live.trees <- nrow(local.live.trees)/nrow(local.trees)
                    
                    condition.hectare <- condition.area/2.471 #convert area from acres to hectares
                    
                    #subset trees into diameter classes for diameter diversity index (regionally specific, classes are identical across forest types)
                    class.1.trees <- local.trees[local.trees$DIA >= 2 & local.trees$DIA < 9.9, ]
                    class.1.trees <- nrow(class.1.trees)/condition.hectare
                    class.2.trees <- local.trees[local.trees$DIA >= 9.9 & local.trees$DIA < 19.7, ]
                    class.2.trees <- nrow(class.2.trees)/condition.hectare
                    class.3.trees <- local.trees[local.trees$DIA >= 19.7 & local.trees$DIA < 39.4, ]
                    class.3.trees <- nrow(class.3.trees)/condition.hectare
                    class.4.trees <- local.trees[local.trees$DIA >= 39.4, ]
                    class.4.trees <- nrow(class.4.trees)/condition.hectare
                    
                    #assign scores to minimum value to avoid errors downstream 
                    tree.score <- 100
                    snag.score <- 100
                    debris.score <- 100
                    diamDiv.score <- 0
                    
                    if(p.live.trees >= 0.1){#begin instructions for if more than 10% of trees are alive
                      
                      if(paz.group %in% "white fir - grand fir"){#begin old growth assessment for (pacific northwest) white/grand fir inside nwfp
                        
                        if(mean(local.live.trees$DIA, na.rm=TRUE) >= .5*29.5){#begin instructions if mean diameter meets minimum threshold for PAZ
                          
                          #assign score for living large trees based on regional regression formula, accounting for breaks in density classes
                          local.live.large.dens <- nrow(local.live.trees[local.live.trees$DIA >= 29.5, ])/condition.hectare #calculate density of live trees above the size threshold
                          if(local.live.large.dens < 14.814){tree.score <- 0 + 3.3751*local.live.large.dens}
                          if(local.live.large.dens >= 14.814 & local.live.large.dens < 41.973){tree.score <- 36.3636 + 0.9205*local.live.large.dens}
                          if(local.live.large.dens >= 41.973){tree.score <- 54.7831 + 0.4816*local.live.large.dens}
                          tree.score <- min(100,tree.score) #cut tree score off at 100, if needed
                          
                          #assign score for snags
                          local.dead.trees <- local.trees[local.trees$STATUSCD == 2, ] #subset to dead trees
                          local.dead.large.dens <- nrow(local.dead.trees[local.dead.trees$DIA >= 19.7, ])/condition.hectare
                          if(local.dead.large.dens < 9.876){snag.score <- 0 + 5.0627*local.live.large.dens}
                          if(local.dead.large.dens >= 9.876 & local.live.large.dens < 19.8088){snag.score <- 25.1429 + 2.5169*local.live.large.dens}
                          if(local.dead.large.dens >= 19.8088){snag.score <- 66.9783 + 0.4049*local.live.large.dens}
                          snag.score <- min(100,snag.score) #cut snag score off at 100, if needed
                          
                          #assign score for downed wood
                          local.woody.debris <- FIA.woodyDebris[FIA.woodyDebris$PLT_CN %in% local.plot$CN & FIA.woodyDebris$CONDID %in% local.condition$CONDID, ] #subset woody debris to local plot and conditoin
                          local.woody.debris <- local.woody.debris[local.woody.debris$TRANSDIA >= 9.8, ] #subset to debris at least 9.8 inches in diameter
                          debris.cover <- sum(local.woody.debris$COVER_PCT, na.rm=TRUE)
                          if(debris.cover < 1.6277){debris.score <- 0 + 30.7181*debris.cover}
                          if(debris.cover >= 1.6277 & debris.cover < 3.4429){debris.score <- 27.5823 + 13.7725*debris.cover}
                          if(debris.cover >= 3.4429){debris.score <- 63.8187 + 3.2476*debris.cover}
                          debris.score <- min(100,debris.score) #cut debris score off at 100, if needed
                          
                          #calculate diameter diversity index
                          class.1.trees <- ifelse(class.1.trees < 200, 0.005*class.1.trees, 1)
                          class.2.trees <- ifelse(class.2.trees < 75, 0.013333*class.2.trees, 1)
                          class.3.trees <- ifelse(class.1.trees < 40, 0.025*class.3.trees, 1)
                          class.4.trees <- ifelse(class.4.trees < 30, 0.03333*class.4.trees, 1)
                          diamDiv.score <- 10*(class.1.trees + (2*class.2.trees) + (3*class.3.trees) + (4*class.4.trees))
                          
                          OGSI200 <- sum(tree.score,snag.score,debris.score,diamDiv.score, na.rm=TRUE)/4
                          
                          local.MOG.status <- ifelse(OGSI200 >= 48, 1, 0) #assign old growth (1) if OFSI200 index surpasses threshold specific to the forest group
                          MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                          
                        }#end instructions if mean diameter meets minimum threshold for PAZ
                      }#end old growth assessment for (pacific northwest) white/grand fir inside nwfp
                      
                      if(paz.group %in% "juniper"){#begin old growth assessment for (pacific northwest) juniper inside nwfp
                        
                        if(mean(local.live.trees$DIA, na.rm=TRUE) >= .5*19.7){#begin instructions if mean diameter meets minimum threshold for PAZ
                          
                          #assign score for living large trees based on regional regression formula, accounting for breaks in density classes
                          local.live.large.dens <- nrow(local.live.trees[local.live.trees$DIA >= 19.7, ])/condition.hectare #calculate density of live trees above the size threshold
                          if(local.live.large.dens < .1){tree.score <- 0 + 500*local.live.large.dens}
                          if(local.live.large.dens >= .1 & local.live.large.dens < 14.8708){tree.score <- 49.8307 + 1.6925*local.live.large.dens}
                          if(local.live.large.dens >= 14.8708){tree.score <- 66.5006 + 0.5715*local.live.large.dens}
                          tree.score <- min(100,tree.score) #cut tree score off at 100, if needed
                          
                          #there is no snag threshold for this forest type
                          #there is no debris threshold for this forest type
                          #there is no diameter diversity threshold for this forest type
                          
                          OGSI200 <- tree.score
                          
                          local.MOG.status <- ifelse(OGSI200 >= 70.95, 1, 0) #assign old growth (1) if OFSI200 index surpasses threshold specific to the forest group
                          MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                          
                        }#end instructions if mean diameter meets minimum threshold for PAZ
                      }#end old growth assessment for (pacific northwest) juniper inside nwfp
                      
                      if(paz.group %in% "mountain hemlock"){#begin old growth assessment for (pacific northwest) mountain hemlock inside nwfp
                        
                        if(mean(local.live.trees$DIA, na.rm=TRUE) >= .5*29.5){#begin instructions if mean diameter meets minimum threshold for PAZ
                          
                          #assign score for living large trees based on regional regression formula, accounting for breaks in density classes
                          local.live.large.dens <- nrow(local.live.trees[local.live.trees$DIA >= 29.5, ])/condition.hectare #calculate density of live trees above the size threshold
                          if(local.live.large.dens < 9.876){tree.score <- 0 + 5.0627*local.live.large.dens}
                          if(local.live.large.dens >= 9.876 & local.live.large.dens < 34.5663){tree.score <- 40.0001 + 1.0125*local.live.large.dens}
                          if(local.live.large.dens >= 34.5663){tree.score <- 56.6049 + 0.5321*local.live.large.dens}
                          tree.score <- min(100,tree.score) #cut tree score off at 100, if needed
                          
                          #assign score for snags
                          local.dead.trees <- local.trees[local.trees$STATUSCD == 2, ] #subset to dead trees
                          local.dead.large.dens <- nrow(local.dead.trees[local.dead.trees$DIA >= 19.7, ])/condition.hectare
                          if(local.dead.large.dens < 12.345){snag.score <- 0 + 4.0502*local.live.large.dens}
                          if(local.dead.large.dens >= 12.345 & local.live.large.dens < 24.7468){snag.score <- 25.1145 + 2.0158*local.live.large.dens}
                          if(local.dead.large.dens >= 24.7468){snag.score <- 62.6278 + 0.4999*local.live.large.dens}
                          snag.score <- min(100,snag.score) #cut snag score off at 100, if needed
                          
                          #assign score for downed wood
                          local.woody.debris <- FIA.woodyDebris[FIA.woodyDebris$PLT_CN %in% local.plot$CN & FIA.woodyDebris$CONDID %in% local.condition$CONDID, ] #subset woody debris to local plot and conditoin
                          local.woody.debris <- local.woody.debris[local.woody.debris$TRANSDIA >= 9.8, ] #subset to debris at least 9.8 inches in diameter
                          debris.cover <- sum(local.woody.debris$COVER_PCT, na.rm=TRUE)
                          if(debris.cover < 2.1817){debris.score <- 0 + 22.9179*debris.cover}
                          if(debris.cover >= 2.1817 & debris.cover < 4.3632){debris.score <- 24.9982 + 11.4597*debris.cover}
                          if(debris.cover >= 4.3632){debris.score <- 59.4239 + 3.5698*debris.cover}
                          debris.score <- min(100,debris.score) #cut debris score off at 100, if needed
                          
                          #calculate diameter diversity index
                          class.1.trees <- ifelse(class.1.trees < 200, 0.005*class.1.trees, 1)
                          class.2.trees <- ifelse(class.2.trees < 75, 0.013333*class.2.trees, 1)
                          class.3.trees <- ifelse(class.1.trees < 40, 0.025*class.3.trees, 1)
                          class.4.trees <- ifelse(class.4.trees < 30, 0.03333*class.4.trees, 1)
                          diamDiv.score <- 10*(class.1.trees + (2*class.2.trees) + (3*class.3.trees) + (4*class.4.trees))
                          
                          OGSI200 <- sum(tree.score,snag.score,debris.score,diamDiv.score, na.rm=TRUE)/4
                          
                          local.MOG.status <- ifelse(OGSI200 >= 41.34, 1, 0) #assign old growth (1) if OFSI200 index surpasses threshold specific to the forest group
                          MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                          
                        }#end instructions if mean diameter meets minimum threshold for PAZ
                      }#end old growth assessment for (pacific northwest) mountain hemlock inside nwfp
                      
                      if(paz.group %in% "oak woodland"){#begin old growth assessment for (pacific northwest) oak/hardwoods inside nwfp
                        
                        if(mean(local.live.trees$DIA, na.rm=TRUE) >= .5*19.7){#begin instructions if mean diameter meets minimum threshold for PAZ
                          
                          #assign score for living large trees based on regional regression formula, accounting for breaks in density classes
                          local.live.large.dens <- nrow(local.live.trees[local.live.trees$DIA >= 19.7, ])/condition.hectare #calculate density of live trees above the size threshold
                          if(local.live.large.dens < 2.469){tree.score <- 0 + 20.2511*local.live.large.dens}
                          if(local.live.large.dens >= 2.469 & local.live.large.dens < 23.0339){tree.score <- 46.9985 + 1.2156*local.live.large.dens}
                          if(local.live.large.dens >= 23.0339){tree.score <- 66.9116 + 0.3511*local.live.large.dens}
                          tree.score <- min(100,tree.score) #cut tree score off at 100, if needed
                          
                          #there is no snag threshold for this forest type
                          #there is no debris threshold for this forest type
                          #there is no diameter diversity threshold for this forest type
                          
                          OGSI200 <- tree.score
                          
                          local.MOG.status <- ifelse(OGSI200 >= 62.46, 1, 0) #assign old growth (1) if OFSI200 index surpasses threshold specific to the forest group
                          MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                          
                        }#end instructions if mean diameter meets minimum threshold for PAZ
                      }#end old growth assessment for (pacific northwest) oak/hardwoods inside nwfp
                      
                      if(paz.group %in% "ponderosa pine"){#begin old growth assessment for (pacific northwest) ponderosa pine forests inside nwfp
                        
                        if(mean(local.live.trees$DIA, na.rm=TRUE) >= .5*29.5){#begin instructions if mean diameter meets minimum threshold for PAZ
                          
                          #assign score for living large trees based on regional regression formula, accounting for breaks in density classes
                          local.live.large.dens <- nrow(local.live.trees[local.live.trees$DIA >= 29.5, ])/condition.hectare #calculate density of live trees above the size threshold
                          if(local.live.large.dens < 2.469){tree.score <- 0 + 20.2511*local.live.large.dens}
                          if(local.live.large.dens >= 2.469 & local.live.large.dens < 12.345){tree.score <- 43.75 + 2.5313*local.live.large.dens}
                          if(local.live.large.dens >= 12.345){tree.score <- 61.1419 + 1.1225*local.live.large.dens}
                          tree.score <- min(100,tree.score) #cut tree score off at 100, if needed
                          
                          #there is no snag threshold for this forest type
                          #there is no debris threshold for this forest type
                          #there is no diameter diversity threshold for this forest type
                          
                          OGSI200 <- tree.score
                          
                          local.MOG.status <- ifelse(OGSI200 >= 67.95, 1, 0) #assign old growth (1) if OFSI200 index surpasses threshold specific to the forest group
                          MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                          
                        }#end instructions if mean diameter meets minimum threshold for PAZ
                      }#end old growth assessment for (pacific northwest) ponderosa pine forests inside nwfp
                      
                      if(paz.group %in% "port orford cedar"){#begin old growth assessment for (pacific northwest) port orford cedar inside nwfp
                        
                        if(mean(local.live.trees$DIA, na.rm=TRUE) >= .5*29.5){#begin instructions if mean diameter meets minimum threshold for PAZ
                          
                          #assign score for living large trees based on regional regression formula, accounting for breaks in density classes
                          local.live.large.dens <- nrow(local.live.trees[local.live.trees$DIA >= 29.5, ])/condition.hectare #calculate density of live trees above the size threshold
                          if(local.live.large.dens < 12.345){tree.score <- 0 + 4.0502*local.live.large.dens}
                          if(local.live.large.dens >= 12.345 & local.live.large.dens < 32.097){tree.score <- 34.375 + 1.2656*local.live.large.dens}
                          if(local.live.large.dens >= 32.097){tree.score <- 58.5609 + 0.5121*local.live.large.dens}
                          tree.score <- min(100,tree.score) #cut tree score off at 100, if needed
                          
                          #assign score for snags
                          local.dead.trees <- local.trees[local.trees$STATUSCD == 2, ] #subset to dead trees
                          local.dead.large.dens <- nrow(local.dead.trees[local.dead.trees$DIA >= 19.7, ])/condition.hectare
                          if(local.dead.large.dens < 14.8708){snag.score <- 0 + 3.3622*local.live.large.dens}
                          if(local.dead.large.dens >= 14.8708 & local.live.large.dens < 31.5365){snag.score <- 27.6925 + 1.5*local.live.large.dens}
                          if(local.dead.large.dens >= 31.5365){snag.score <- 30.0436 + 1.4255*local.live.large.dens}
                          snag.score <- min(100,snag.score) #cut snag score off at 100, if needed
                          
                          #assign score for downed wood
                          local.woody.debris <- FIA.woodyDebris[FIA.woodyDebris$PLT_CN %in% local.plot$CN & FIA.woodyDebris$CONDID %in% local.condition$CONDID, ] #subset woody debris to local plot and conditoin
                          local.woody.debris <- local.woody.debris[local.woody.debris$TRANSDIA >= 9.8, ] #subset to debris at least 9.8 inches in diameter
                          debris.cover <- sum(local.woody.debris$COVER_PCT, na.rm=TRUE)
                          if(debris.cover < 0.9715){debris.score <- 0 + 51.4641*debris.cover}
                          if(debris.cover >= 0.9715 & debris.cover < 2.6162){debris.score <- 35.2318 + 15.2005*debris.cover}
                          if(debris.cover >= 2.6162){debris.score <- 61.9655 + 4.9821*debris.cover}
                          debris.score <- min(100,debris.score) #cut debris score off at 100, if needed
                          
                          #calculate diameter diversity index
                          class.1.trees <- ifelse(class.1.trees < 200, 0.005*class.1.trees, 1)
                          class.2.trees <- ifelse(class.2.trees < 75, 0.013333*class.2.trees, 1)
                          class.3.trees <- ifelse(class.1.trees < 40, 0.025*class.3.trees, 1)
                          class.4.trees <- ifelse(class.4.trees < 30, 0.03333*class.4.trees, 1)
                          diamDiv.score <- 10*(class.1.trees + (2*class.2.trees) + (3*class.3.trees) + (4*class.4.trees))
                          
                          OGSI200 <- sum(tree.score,snag.score,debris.score,diamDiv.score, na.rm=TRUE)/4
                          
                          local.MOG.status <- ifelse(OGSI200 >= 45.01, 1, 0) #assign old growth (1) if OFSI200 index surpasses threshold specific to the forest group
                          MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                          
                        }#end instructions if mean diameter meets minimum threshold for PAZ
                      }#end old growth assessment for (pacific northwest) port orford cedar inside nwfp
                      
                      if(paz.group %in% "shasta red fir"){#begin old growth assessment for (pacific northwest) California shasta and red fir inside nwfp
                        
                        if(mean(local.live.trees$DIA, na.rm=TRUE) >= .5*29.5){#begin instructions if mean diameter meets minimum threshold for PAZ
                          
                          #assign score for living large trees based on regional regression formula, accounting for breaks in density classes
                          local.live.large.dens <- nrow(local.live.trees[local.live.trees$DIA >= 29.5, ])/condition.hectare #calculate density of live trees above the size threshold
                          if(local.live.large.dens < 19.752){tree.score <- 0 + 2.5313*local.live.large.dens}
                          if(local.live.large.dens >= 19.752 & local.live.large.dens < 46.911){tree.score <- 31.8181 + 0.9205*local.live.large.dens}
                          if(local.live.large.dens >= 46.911){tree.score <- 38.4615 + 0.7788*local.live.large.dens}
                          tree.score <- min(100,tree.score) #cut tree score off at 100, if needed
                          
                          #assign score for snags
                          local.dead.trees <- local.trees[local.trees$STATUSCD == 2, ] #subset to dead trees
                          local.dead.large.dens <- nrow(local.dead.trees[local.dead.trees$DIA >= 19.7, ])/condition.hectare
                          if(local.dead.large.dens < 7.407){snag.score <- 0 + 6.7503*local.live.large.dens}
                          if(local.dead.large.dens >= 7.407 & local.live.large.dens < 14.814){snag.score <- 25 + 3.3751*local.live.large.dens}
                          if(local.dead.large.dens >= 14.814){snag.score <- 64.8879 + 0.6826*local.live.large.dens}
                          snag.score <- min(100,snag.score) #cut snag score off at 100, if needed
                          
                          #assign score for downed wood
                          local.woody.debris <- FIA.woodyDebris[FIA.woodyDebris$PLT_CN %in% local.plot$CN & FIA.woodyDebris$CONDID %in% local.condition$CONDID, ] #subset woody debris to local plot and conditoin
                          local.woody.debris <- local.woody.debris[local.woody.debris$TRANSDIA >= 9.8, ] #subset to debris at least 9.8 inches in diameter
                          debris.cover <- sum(local.woody.debris$COVER_PCT, na.rm=TRUE)
                          if(debris.cover < 0.6477){debris.score <- 0 + 77.1962*debris.cover}
                          if(debris.cover >= 0.6477 & debris.cover < 2.2328){debris.score <- 39.7847 + 15.7716*debris.cover}
                          if(debris.cover >= 2.2328){debris.score <- 66.1432 + 3.9666*debris.cover}
                          debris.score <- min(100,debris.score) #cut debris score off at 100, if needed
                          
                          #calculate diameter diversity index
                          class.1.trees <- ifelse(class.1.trees < 200, 0.005*class.1.trees, 1)
                          class.2.trees <- ifelse(class.2.trees < 75, 0.013333*class.2.trees, 1)
                          class.3.trees <- ifelse(class.1.trees < 40, 0.025*class.3.trees, 1)
                          class.4.trees <- ifelse(class.4.trees < 30, 0.03333*class.4.trees, 1)
                          diamDiv.score <- 10*(class.1.trees + (2*class.2.trees) + (3*class.3.trees) + (4*class.4.trees))
                          
                          OGSI200 <- sum(tree.score,snag.score,debris.score,diamDiv.score, na.rm=TRUE)/4
                          
                          local.MOG.status <- ifelse(OGSI200 >= 52.81, 1, 0) #assign old growth (1) if OFSI200 index surpasses threshold specific to the forest group
                          MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                          
                        }#end instructions if mean diameter meets minimum threshold for PAZ
                      }#end old growth assessment for (pacific northwest) California shasta and red fir inside nwfp
                      
                      if(paz.group %in% "silver fir"){#begin old growth assessment for (pacific northwest) silver fir inside nwfp
                        
                        if(mean(local.live.trees$DIA, na.rm=TRUE) >= .5*29.5){#begin instructions if mean diameter meets minimum threshold for PAZ
                          
                          #assign score for living large trees based on regional regression formula, accounting for breaks in density classes
                          local.live.large.dens <- nrow(local.live.trees[local.live.trees$DIA >= 29.5, ])/condition.hectare #calculate density of live trees above the size threshold
                          if(local.live.large.dens < 24.2677){tree.score <- 0 + 2.0603*local.live.large.dens}
                          if(local.live.large.dens >= 24.2677 & local.live.large.dens < 54.318){tree.score <- 29.8106 + 0.8319*local.live.large.dens}
                          if(local.live.large.dens >= 54.318){tree.score <- 46.3446 + 0.5275*local.live.large.dens}
                          tree.score <- min(100,tree.score) #cut tree score off at 100, if needed
                          
                          #assign score for snags
                          local.dead.trees <- local.trees[local.trees$STATUSCD == 2, ] #subset to dead trees
                          local.dead.large.dens <- nrow(local.dead.trees[local.dead.trees$DIA >= 19.7, ])/condition.hectare
                          if(local.dead.large.dens < 21.42){snag.score <- 0 + 2.3342*local.live.large.dens}
                          if(local.dead.large.dens >= 21.42 & local.live.large.dens < 37.1711){snag.score <- 16.0024 + 1.5871*local.live.large.dens}
                          if(local.dead.large.dens >= 37.1711){snag.score <- 55.5568 + 0.523*local.live.large.dens}
                          snag.score <- min(100,snag.score) #cut snag score off at 100, if needed
                          
                          #assign score for downed wood
                          local.woody.debris <- FIA.woodyDebris[FIA.woodyDebris$PLT_CN %in% local.plot$CN & FIA.woodyDebris$CONDID %in% local.condition$CONDID, ] #subset woody debris to local plot and conditoin
                          local.woody.debris <- local.woody.debris[local.woody.debris$TRANSDIA >= 9.8, ] #subset to debris at least 9.8 inches in diameter
                          debris.cover <- sum(local.woody.debris$COVER_PCT, na.rm=TRUE)
                          if(debris.cover < 4.3974){debris.score <- 0 + 11.3703*debris.cover}
                          if(debris.cover >= 4.3974 & debris.cover < 7.7551){debris.score <- 17.2593 + 7.4454*debris.cover}
                          if(debris.cover >= 7.7551){debris.score <- 60.9019 + 1.8179*debris.cover}
                          debris.score <- min(100,debris.score) #cut debris score off at 100, if needed
                          
                          #calculate diameter diversity index
                          class.1.trees <- ifelse(class.1.trees < 200, 0.005*class.1.trees, 1)
                          class.2.trees <- ifelse(class.2.trees < 75, 0.013333*class.2.trees, 1)
                          class.3.trees <- ifelse(class.1.trees < 40, 0.025*class.3.trees, 1)
                          class.4.trees <- ifelse(class.4.trees < 30, 0.03333*class.4.trees, 1)
                          diamDiv.score <- 10*(class.1.trees + (2*class.2.trees) + (3*class.3.trees) + (4*class.4.trees))
                          
                          OGSI200 <- sum(tree.score,snag.score,debris.score,diamDiv.score, na.rm=TRUE)/4
                          
                          local.MOG.status <- ifelse(OGSI200 >= 43.39, 1, 0) #assign old growth (1) if OFSI200 index surpasses threshold specific to the forest group
                          MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                          
                        }#end instructions if mean diameter meets minimum threshold for PAZ
                      }#end old growth assessment for (pacific northwest) silver fir inside nwfp
                      
                      if(paz.group %in% "sitka spruce"){#begin old growth assessment for (pacific northwest) sitka spruce inside nwfp
                        
                        if(mean(local.live.trees$DIA, na.rm=TRUE) >= .5*39.4){#begin instructions if mean diameter meets minimum threshold for PAZ
                          
                          #assign score for living large trees based on regional regression formula, accounting for breaks in density classes
                          local.live.large.dens <- nrow(local.live.trees[local.live.trees$DIA >= 39.4, ])/condition.hectare #calculate density of live trees above the size threshold
                          if(local.live.large.dens < 9.2587){tree.score <- 0 + 5.4002*local.live.large.dens}
                          if(local.live.large.dens >= 9.2587 & local.live.large.dens < 28.3935){tree.score <- 37.9032 + 1.3065*local.live.large.dens}
                          if(local.live.large.dens >= 28.3935){tree.score <- 44.086 + 1.0887*local.live.large.dens}
                          tree.score <- min(100,tree.score) #cut tree score off at 100, if needed
                          
                          #assign score for snags
                          local.dead.trees <- local.trees[local.trees$STATUSCD == 2, ] #subset to dead trees
                          local.dead.large.dens <- nrow(local.dead.trees[local.dead.trees$DIA >= 19.7, ])/condition.hectare
                          if(local.dead.large.dens < 9.876){snag.score <- 0 + 5.0627*local.live.large.dens}
                          if(local.dead.large.dens >= 9.876 & local.live.large.dens < 16.0485){snag.score <- 10 + 4.0502*local.live.large.dens}
                          if(local.dead.large.dens >= 16.0485){snag.score <- 51.4492 + 1.4674*local.live.large.dens}
                          snag.score <- min(100,snag.score) #cut snag score off at 100, if needed
                          
                          #assign score for downed wood
                          local.woody.debris <- FIA.woodyDebris[FIA.woodyDebris$PLT_CN %in% local.plot$CN & FIA.woodyDebris$CONDID %in% local.condition$CONDID, ] #subset woody debris to local plot and conditoin
                          local.woody.debris <- local.woody.debris[local.woody.debris$TRANSDIA >= 9.8, ] #subset to debris at least 9.8 inches in diameter
                          debris.cover <- sum(local.woody.debris$COVER_PCT, na.rm=TRUE)
                          if(debris.cover < 4.6871){debris.score <- 0 + 10.6674*debris.cover}
                          if(debris.cover >= 4.6871 & debris.cover < 7.1243){debris.score <- 1.9207 + 10.2576*debris.cover}
                          if(debris.cover >= 7.1243){debris.score <- 53.8812 + 2.9643*debris.cover}
                          debris.score <- min(100,debris.score) #cut debris score off at 100, if needed
                          
                          #calculate diameter diversity index
                          class.1.trees <- ifelse(class.1.trees < 200, 0.005*class.1.trees, 1)
                          class.2.trees <- ifelse(class.2.trees < 75, 0.013333*class.2.trees, 1)
                          class.3.trees <- ifelse(class.1.trees < 40, 0.025*class.3.trees, 1)
                          class.4.trees <- ifelse(class.4.trees < 30, 0.03333*class.4.trees, 1)
                          diamDiv.score <- 10*(class.1.trees + (2*class.2.trees) + (3*class.3.trees) + (4*class.4.trees))
                          
                          OGSI200 <- sum(tree.score,snag.score,debris.score,diamDiv.score, na.rm=TRUE)/4
                          
                          local.MOG.status <- ifelse(OGSI200 >= 59.96, 1, 0) #assign old growth (1) if OFSI200 index surpasses threshold specific to the forest group
                          MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                          
                        }#end instructions if mean diameter meets minimum threshold for PAZ
                      }#end old growth assessment for (pacific northwest) sitka spruce inside nwfp
                      
                      if(paz.group %in% "subalpine fir"){#begin old growth assessment for (pacific northwest) subalpine fir and engelmann spruce inside nwfp
                        
                        if(mean(local.live.trees$DIA, na.rm=TRUE) >= .5*19.7){#begin instructions if mean diameter meets minimum threshold for PAZ
                          
                          #assign score for living large trees based on regional regression formula, accounting for breaks in density classes
                          local.live.large.dens <- nrow(local.live.trees[local.live.trees$DIA >= 19.7, ])/condition.hectare #calculate density of live trees above the size threshold
                          if(local.live.large.dens < 14.814){tree.score <- 0 + 3.3751*local.live.large.dens}
                          if(local.live.large.dens >= 14.814 & local.live.large.dens < 59.3242){tree.score <- 41.6794 + 0.5616*local.live.large.dens}
                          if(local.live.large.dens >= 59.3242){tree.score <- 55.1433 + 0.3347*local.live.large.dens}
                          tree.score <- min(100,tree.score) #cut tree score off at 100, if needed
                          
                          #assign score for snags
                          local.dead.trees <- local.trees[local.trees$STATUSCD == 2, ] #subset to dead trees
                          local.dead.large.dens <- nrow(local.dead.trees[local.dead.trees$DIA >= 19.7, ])/condition.hectare
                          if(local.dead.large.dens < 2.6195){snag.score <- 0 + 19.0876*local.live.large.dens}
                          if(local.dead.large.dens >= 2.6195 & local.live.large.dens < 14.8708){snag.score <- 44.6546 + 2.0406*local.live.large.dens}
                          if(local.dead.large.dens >= 14.8708){snag.score <- 68.6455 + 0.4273*local.live.large.dens}
                          snag.score <- min(100,snag.score) #cut snag score off at 100, if needed
                          
                          #assign score for downed wood
                          local.woody.debris <- FIA.woodyDebris[FIA.woodyDebris$PLT_CN %in% local.plot$CN & FIA.woodyDebris$CONDID %in% local.condition$CONDID, ] #subset woody debris to local plot and conditoin
                          local.woody.debris <- local.woody.debris[local.woody.debris$TRANSDIA >= 9.8, ] #subset to debris at least 9.8 inches in diameter
                          debris.cover <- sum(local.woody.debris$COVER_PCT, na.rm=TRUE)
                          if(debris.cover < 1.2612){debris.score <- 0 + 39.6447*debris.cover}
                          if(debris.cover >= 1.2612 & debris.cover < 3.0509){debris.score <- 32.383 + 13.9684*debris.cover}
                          if(debris.cover >= 3.0509){debris.score <- 66.5385 + 2.7733*debris.cover}
                          debris.score <- min(100,debris.score) #cut debris score off at 100, if needed
                          
                          #calculate diameter diversity index
                          class.1.trees <- ifelse(class.1.trees < 200, 0.005*class.1.trees, 1)
                          class.2.trees <- ifelse(class.2.trees < 75, 0.013333*class.2.trees, 1)
                          class.3.trees <- ifelse(class.1.trees < 40, 0.025*class.3.trees, 1)
                          class.4.trees <- ifelse(class.4.trees < 30, 0.03333*class.4.trees, 1)
                          diamDiv.score <- 10*(class.1.trees + (2*class.2.trees) + (3*class.3.trees) + (4*class.4.trees))
                          
                          OGSI200 <- sum(tree.score,snag.score,debris.score,diamDiv.score, na.rm=TRUE)/4
                          
                          local.MOG.status <- ifelse(OGSI200 >= 45.08, 1, 0) #assign old growth (1) if OFSI200 index surpasses threshold specific to the forest group
                          MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                          
                        }#end instructions if mean diameter meets minimum threshold for PAZ
                      }#end old growth assessment for (pacific northwest) subalpine fir and engelmann spruce inside nwfp
                      
                      if(paz.group %in% "tanoak"){#begin old growth assessment for (pacific northwest) tanoak inside nwfp
                        
                        if(mean(local.live.trees$DIA, na.rm=TRUE) >= .5*39.4){#begin instructions if mean diameter meets minimum threshold for PAZ
                          
                          #assign score for living large trees based on regional regression formula, accounting for breaks in density classes
                          local.live.large.dens <- nrow(local.live.trees[local.live.trees$DIA >= 39.4, ])/condition.hectare #calculate density of live trees above the size threshold
                          if(local.live.large.dens < 12.345){tree.score <- 0 + 4.0502*local.live.large.dens}
                          if(local.live.large.dens >= 12.345 & local.live.large.dens < 24.69){tree.score <- 25 + 2.0251*local.live.large.dens}
                          if(local.live.large.dens >= 24.69){tree.score <- 46.7832 + 1.1428*local.live.large.dens}
                          tree.score <- min(100,tree.score) #cut tree score off at 100, if needed
                          
                          #assign score for snags
                          local.dead.trees <- local.trees[local.trees$STATUSCD == 2, ] #subset to dead trees
                          local.dead.large.dens <- nrow(local.dead.trees[local.dead.trees$DIA >= 39.4, ])/condition.hectare
                          if(local.dead.large.dens < 4.938){snag.score <- 0 + 10.1255*local.live.large.dens}
                          if(local.dead.large.dens >= 4.938 & local.live.large.dens < 14.9374){snag.score <- 37.6542 + 2.5001*local.live.large.dens}
                          if(local.dead.large.dens >= 14.9374){snag.score <- 64.8584 + 0.6789*local.live.large.dens}
                          snag.score <- min(100,snag.score) #cut snag score off at 100, if needed
                          
                          #assign score for downed wood
                          local.woody.debris <- FIA.woodyDebris[FIA.woodyDebris$PLT_CN %in% local.plot$CN & FIA.woodyDebris$CONDID %in% local.condition$CONDID, ] #subset woody debris to local plot and conditoin
                          local.woody.debris <- local.woody.debris[local.woody.debris$TRANSDIA >= 9.8, ] #subset to debris at least 9.8 inches in diameter
                          debris.cover <- sum(local.woody.debris$COVER_PCT, na.rm=TRUE)
                          if(debris.cover < 1.3635){debris.score <- 0 + 36.6703*debris.cover}
                          if(debris.cover >= 1.3635 & debris.cover < 3.4982){debris.score <- 34.032 + 11.7109*debris.cover}
                          if(debris.cover >= 3.4982){debris.score <- 61.7148 + 3.7976*debris.cover}
                          debris.score <- min(100,debris.score) #cut debris score off at 100, if needed
                          
                          #calculate diameter diversity index
                          class.1.trees <- ifelse(class.1.trees < 200, 0.005*class.1.trees, 1)
                          class.2.trees <- ifelse(class.2.trees < 75, 0.013333*class.2.trees, 1)
                          class.3.trees <- ifelse(class.1.trees < 40, 0.025*class.3.trees, 1)
                          class.4.trees <- ifelse(class.4.trees < 30, 0.03333*class.4.trees, 1)
                          diamDiv.score <- 10*(class.1.trees + (2*class.2.trees) + (3*class.3.trees) + (4*class.4.trees))
                          
                          OGSI200 <- sum(tree.score,snag.score,debris.score,diamDiv.score, na.rm=TRUE)/4
                          
                          local.MOG.status <- ifelse(OGSI200 >= 48.22, 1, 0) #assign old growth (1) if OFSI200 index surpasses threshold specific to the forest group
                          MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                          
                        }#end instructions if mean diameter meets minimum threshold for PAZ
                      }#end old growth assessment for (pacific northwest) tanoak inside nwfp
                      
                      if(paz.group %in% "western hemlock"){#begin old growth assessment for (pacific northwest) western hemlock inside nwfp
                        
                        if(mean(local.live.trees$DIA, na.rm=TRUE) >= .5*39.4){#begin instructions if mean diameter meets minimum threshold for PAZ
                          
                          #assign score for living large trees based on regional regression formula, accounting for breaks in density classes
                          local.live.large.dens <- nrow(local.live.trees[local.live.trees$DIA >= 39.4, ])/condition.hectare #calculate density of live trees above the size threshold
                          if(local.live.large.dens < 9.876){tree.score <- 0 + 5.0627*local.live.large.dens}
                          if(local.live.large.dens >= 9.876 & local.live.large.dens < 27.254){tree.score <- 35.7923 + 1.4386*local.live.large.dens}
                          if(local.live.large.dens >= 27.254){tree.score <- 55.2341 + 0.7252*local.live.large.dens}
                          tree.score <- min(100,tree.score) #cut tree score off at 100, if needed
                          
                          #assign score for snags
                          local.dead.trees <- local.trees[local.trees$STATUSCD == 2, ] #subset to dead trees
                          local.dead.large.dens <- nrow(local.dead.trees[local.dead.trees$DIA >= 39.4, ])/condition.hectare
                          if(local.dead.large.dens < 6.7116){snag.score <- 0 + 7.4497*local.live.large.dens}
                          if(local.dead.large.dens >= 6.7116 & local.live.large.dens < 12.345){snag.score <- 20.2151 + 4.4378*local.live.large.dens}
                          if(local.dead.large.dens >= 12.345){snag.score <- 62.5286 + 1.0102*local.live.large.dens}
                          snag.score <- min(100,snag.score) #cut snag score off at 100, if needed
                          
                          #assign score for downed wood
                          local.woody.debris <- FIA.woodyDebris[FIA.woodyDebris$PLT_CN %in% local.plot$CN & FIA.woodyDebris$CONDID %in% local.condition$CONDID, ] #subset woody debris to local plot and conditoin
                          local.woody.debris <- local.woody.debris[local.woody.debris$TRANSDIA >= 9.8, ] #subset to debris at least 9.8 inches in diameter
                          debris.cover <- sum(local.woody.debris$COVER_PCT, na.rm=TRUE)
                          if(debris.cover < 3.9884){debris.score <- 0 + 12.5363*debris.cover}
                          if(debris.cover >= 3.9884 & debris.cover < 7.1244){debris.score <- 18.2047 + 7.9719*debris.cover}
                          if(debris.cover >= 7.1244){debris.score <- 59.9428 + 2.1134*debris.cover}
                          debris.score <- min(100,debris.score) #cut debris score off at 100, if needed
                          
                          #calculate diameter diversity index
                          class.1.trees <- ifelse(class.1.trees < 200, 0.005*class.1.trees, 1)
                          class.2.trees <- ifelse(class.2.trees < 75, 0.013333*class.2.trees, 1)
                          class.3.trees <- ifelse(class.1.trees < 40, 0.025*class.3.trees, 1)
                          class.4.trees <- ifelse(class.4.trees < 30, 0.03333*class.4.trees, 1)
                          diamDiv.score <- 10*(class.1.trees + (2*class.2.trees) + (3*class.3.trees) + (4*class.4.trees))
                          
                          OGSI200 <- sum(tree.score,snag.score,debris.score,diamDiv.score, na.rm=TRUE)/4
                          
                          local.MOG.status <- ifelse(OGSI200 >= 44.63, 1, 0) #assign old growth (1) if OFSI200 index surpasses threshold specific to the forest group
                          MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                          
                        }#end instructions if mean diameter meets minimum threshold for PAZ
                      }#end old growth assessment for (pacific northwest) western hemlock inside nwfp
                      
                      if(paz.group %in% "douglas fir"){#begin old growth assessment for (pacific northwest) douglas fir inside nwfp
                        
                        if(mean(local.live.trees$DIA, na.rm=TRUE) >= .5*29.5){#begin instructions if mean diameter meets minimum threshold for PAZ
                          
                          #assign score for living large trees based on regional regression formula, accounting for breaks in density classes
                          local.live.large.dens <- nrow(local.live.trees[local.live.trees$DIA >= 29.5, ])/condition.hectare #calculate density of live trees above the size threshold
                          if(local.live.large.dens < 4.938){tree.score <- 0 + 10.1255*local.live.large.dens}
                          if(local.live.large.dens >= 4.938 & local.live.large.dens < 24.69){tree.score <- 43.75 + 1.2656*local.live.large.dens}
                          if(local.live.large.dens >= 24.69){tree.score <- 59.5297 + 0.6265*local.live.large.dens}
                          tree.score <- min(100,tree.score) #cut tree score off at 100, if needed
                          
                          #assign score for snags
                          local.dead.trees <- local.trees[local.trees$STATUSCD == 2, ] #subset to dead trees
                          local.dead.large.dens <- nrow(local.dead.trees[local.dead.trees$DIA >= 19.7, ])/condition.hectare
                          if(local.dead.large.dens < 2.469){snag.score <- 0 + 20.2511*local.live.large.dens}
                          if(local.dead.large.dens >= 2.469 & local.live.large.dens < 7.5099){snag.score <- 37.7551 + 4.9594*local.live.large.dens}
                          if(local.dead.large.dens >= 7.5099){snag.score <- 70.3665 + 0.6169*local.live.large.dens}
                          snag.score <- min(100,snag.score) #cut snag score off at 100, if needed
                          
                          #assign score for downed wood
                          local.woody.debris <- FIA.woodyDebris[FIA.woodyDebris$PLT_CN %in% local.plot$CN & FIA.woodyDebris$CONDID %in% local.condition$CONDID, ] #subset woody debris to local plot and conditoin
                          local.woody.debris <- local.woody.debris[local.woody.debris$TRANSDIA >= 9.8, ] #subset to debris at least 9.8 inches in diameter
                          debris.cover <- sum(local.woody.debris$COVER_PCT, na.rm=TRUE)
                          if(debris.cover < 0.784){debris.score <- 0 + 63.7755*debris.cover}
                          if(debris.cover >= 0.784 & debris.cover < 1.909){debris.score <- 32.5777 + 22.2222*debris.cover}
                          if(debris.cover >= 1.909){debris.score <- 65.8906 + 4.7717*debris.cover}
                          debris.score <- min(100,debris.score) #cut debris score off at 100, if needed
                          
                          #calculate diameter diversity index
                          class.1.trees <- ifelse(class.1.trees < 200, 0.005*class.1.trees, 1)
                          class.2.trees <- ifelse(class.2.trees < 75, 0.013333*class.2.trees, 1)
                          class.3.trees <- ifelse(class.1.trees < 40, 0.025*class.3.trees, 1)
                          class.4.trees <- ifelse(class.4.trees < 30, 0.03333*class.4.trees, 1)
                          diamDiv.score <- 10*(class.1.trees + (2*class.2.trees) + (3*class.3.trees) + (4*class.4.trees))
                          
                          OGSI200 <- sum(tree.score,snag.score,debris.score,diamDiv.score, na.rm=TRUE)/4
                          
                          local.MOG.status <- ifelse(OGSI200 >= 50.81, 1, 0) #assign old growth (1) if OFSI200 index surpasses threshold specific to the forest group
                          MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                          
                        }#end instructions if mean diameter meets minimum threshold for PAZ
                      }#end old growth assessment for (pacific northwest) douglas fir inside nwfp
                      
                      if(paz.group %in% "lodgepole pine"){#begin old growth assessment for (pacific northwest) lodgepole pine forests inside nwfp
                        
                        if(mean(local.live.trees$DIA, na.rm=TRUE) >= .5*9.8){#begin instructions if mean diameter meets minimum threshold for PAZ
                          
                          #assign score for living large trees based on regional regression formula, accounting for breaks in density classes
                          local.live.large.dens <- nrow(local.live.trees[local.live.trees$DIA >= 9.8, ])/condition.hectare #calculate density of live trees above the size threshold
                          if(local.live.large.dens < 55.984){tree.score <- 0 + 0.8931*local.live.large.dens}
                          if(local.live.large.dens >= 55.984 & local.live.large.dens < 247.752){tree.score <- 42.7016 + 0.1303*local.live.large.dens}
                          if(local.live.large.dens >= 247.752){tree.score <- 51.7325 + 0.0939*local.live.large.dens}
                          tree.score <- min(100,tree.score) #cut tree score off at 100, if needed
                          
                          #there is no snag threshold for this forest type
                          #there is no debris threshold for this forest type
                          #there is no diameter diversity threshold for this forest type
                          
                          OGSI200 <- tree.score
                          
                          local.MOG.status <- ifelse(OGSI200 >= 61.46, 1, 0) #assign old growth (1) if OFSI200 index surpasses threshold specific to the forest group
                          MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                          
                        }#end instructions if mean diameter meets minimum threshold for PAZ
                      }#end old growth assessment for (pacific northwest) lodgepole pine forests inside nwfp
                      
                      if(paz.group %in% "jeffrey pine - knobcone pine"){#begin old growth assessment for (pacific northwest) jeffrey/knobcone pine forests inside nwfp
                        
                        if(mean(local.live.trees$DIA, na.rm=TRUE) >= .5*29.5){#begin instructions if mean diameter meets minimum threshold for PAZ
                          
                          #assign score for living large trees based on regional regression formula, accounting for breaks in density classes
                          local.live.large.dens <- nrow(local.live.trees[local.live.trees$DIA >= 29.5, ])/condition.hectare #calculate density of live trees above the size threshold
                          if(local.live.large.dens < 4.938){tree.score <- 0 + 10.1255*local.live.large.dens}
                          if(local.live.large.dens >= 4.938 & local.live.large.dens < 17.283){tree.score <- 40 + 2.0251*local.live.large.dens}
                          if(local.live.large.dens >= 17.283){tree.score <- 57.9434 + 0.9868*local.live.large.dens}
                          tree.score <- min(100,tree.score) #cut tree score off at 100, if needed
                          
                          #there is no snag threshold for this forest type
                          #there is no debris threshold for this forest type
                          #there is no diameter diversity threshold for this forest type
                          
                          OGSI200 <- tree.score
                          
                          local.MOG.status <- ifelse(OGSI200 >= 61.77, 1, 0) #assign old growth (1) if OFSI200 index surpasses threshold specific to the forest group
                          MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                          
                        }#end instructions if mean diameter meets minimum threshold for PAZ
                      }#end old growth assessment for (pacific northwest) jeffrey/knobcone pine forests inside nwfp
                      
                      if(paz.group %in% "redwood"){#begin old growth assessment for (pacific northwest) redwoods inside nwfp
                        
                        if(mean(local.live.trees$DIA, na.rm=TRUE) >= .5*39.4){#begin instructions if mean diameter meets minimum threshold for PAZ
                          
                          #assign score for living large trees based on regional regression formula, accounting for breaks in density classes
                          local.live.large.dens <- nrow(local.live.trees[local.live.trees$DIA >= 39.4, ])/condition.hectare #calculate density of live trees above the size threshold
                          if(local.live.large.dens < 16.6657){tree.score <- 0 + 3.0001*local.live.large.dens}
                          if(local.live.large.dens >= 16.6657 & local.live.large.dens < 37.6522){tree.score <- 30.147 + 1.1912*local.live.large.dens}
                          if(local.live.large.dens >= 37.6522){tree.score <- 47.2727 + 0.7364*local.live.large.dens}
                          tree.score <- min(100,tree.score) #cut tree score off at 100, if needed
                          
                          #assign score for snags
                          local.dead.trees <- local.trees[local.trees$STATUSCD == 2, ] #subset to dead trees
                          local.dead.large.dens <- nrow(local.dead.trees[local.dead.trees$DIA >= 39.4, ])/condition.hectare
                          if(local.dead.large.dens < 2.469){snag.score <- 0 + 20.2511*local.live.large.dens}
                          if(local.dead.large.dens >= 2.469 & local.live.large.dens < 5.5319){snag.score <- 29.8478 + 8.162*local.live.large.dens}
                          if(local.dead.large.dens >= 5.5319){snag.score <- 63.2309 + 2.1274*local.live.large.dens}
                          snag.score <- min(100,snag.score) #cut snag score off at 100, if needed
                          
                          #assign score for downed wood
                          local.woody.debris <- FIA.woodyDebris[FIA.woodyDebris$PLT_CN %in% local.plot$CN & FIA.woodyDebris$CONDID %in% local.condition$CONDID, ] #subset woody debris to local plot and conditoin
                          local.woody.debris <- local.woody.debris[local.woody.debris$TRANSDIA >= 9.8, ] #subset to debris at least 9.8 inches in diameter
                          debris.cover <- sum(local.woody.debris$COVER_PCT, na.rm=TRUE)
                          if(debris.cover < 2.0109){debris.score <- 0 + 24.8644*debris.cover}
                          if(debris.cover >= 2.0109 & debris.cover < 4.5045){debris.score <- 29.8397 + 10.0254*debris.cover}
                          if(debris.cover >= 4.5045){debris.score <- 65.7197 + 2.0601*debris.cover}
                          debris.score <- min(100,debris.score) #cut debris score off at 100, if needed
                          
                          #calculate diameter diversity index
                          class.1.trees <- ifelse(class.1.trees < 200, 0.005*class.1.trees, 1)
                          class.2.trees <- ifelse(class.2.trees < 75, 0.013333*class.2.trees, 1)
                          class.3.trees <- ifelse(class.1.trees < 40, 0.025*class.3.trees, 1)
                          class.4.trees <- ifelse(class.4.trees < 30, 0.03333*class.4.trees, 1)
                          diamDiv.score <- 10*(class.1.trees + (2*class.2.trees) + (3*class.3.trees) + (4*class.4.trees))
                          
                          OGSI200 <- sum(tree.score,snag.score,debris.score,diamDiv.score, na.rm=TRUE)/4
                          
                          local.MOG.status <- ifelse(OGSI200 >= 50.94, 1, 0) #assign old growth (1) if OFSI200 index surpasses threshold specific to the forest group
                          MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                          
                        }#end instructions if mean diameter meets minimum threshold for PAZ
                      }#end old growth assessment for (pacific northwest) redwoods inside nwfp
                      
                    }#end instructions for if more than 10% of trees are alive
                  }#end instructions for if the plot is inside the northwest forest plan
                  
                  #begin (pacific northwest) old growth assessment for areas OUTSIDE of the northwest forest per Table 14 of the MOG document.
                  if(nrow(nwfp.test) == 0){#begin instructions for if the plot is outside the northwest forest plan
                    
                    site.class <- max(local.condition$SITECLCD,local.condition$SITECLCDEST,na.rm=TRUE) #use either the measured or estimated value of productivity (one of these should always be NA, or else they should be the same value)
                    
                    if(paz.group %in% "white fir - grand fir"){#begin old growth instructions for (pacific northwest) white fir and grand fir forests outside of the management zone
                      #Non-NWFP White/grand fir is the only group that straddles two different geographic ranges. Determine which range this is
                      ECOSUBCD <- FIA.geomplot[which(FIA.geomplot$CN %in% local.plot$CN),"ECOSUBCD"] #assign ecological subregion
                      regional.geog <- ifelse(substr(ECOSUBCD,1,5) %in% "M242C", "central", "non-central")
                      county.test <- st_filter(x = local.plot, 
                                               y = counties)
                      if(nrow(county.test)>0 & regional.geog %in% "central"){regional.geog <- "non-central"} #if plot falls in Hood River or Wasco Counties, it is not considered central Oregon regardless of its ecosystem classification
                      
                      if(regional.geog %in% "central"){#begin sub-assessment for white fir/grand fir in central part of region
                        
                        if(site.class %in% c(1,2)){#begin sub-sub-assessment for if site class is high
                          
                          local.large.trees <- local.trees[local.trees$DIA >= 21, ] #isolate to large trees specific to regional and forest type threshold
                          large.trees.p.acre <- nrow(local.large.trees)/condition.area #convert from count to n per acre
                          large.trees.p.acre <- ifelse(large.trees.p.acre >= 15, 1, 0)
                          
                          max.age <- max(c(raw.stand.age,local.trees$tree.age), na.rm=TRUE) #determine the age of the oldest tree in the stand
                          max.age <- ifelse(max.age >= 150, 1, 0)
                          
                          local.MOG.status <- ifelse(sum(large.trees.p.acre, max.age, na.rm=TRUE) == 2, 1, 0) #assign old growth status (1) if threshold is met
                          MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                          
                        }#end sub-sub-assessment for if site class is high
                        
                        if(!(site.class %in% c(1,2))){#begin sub-sub-assessment for if site class is low or medium
                          
                          local.large.trees <- local.trees[local.trees$DIA >= 21, ] #isolate to large trees specific to regional and forest type threshold
                          large.trees.p.acre <- nrow(local.large.trees)/condition.area #convert from count to n per acre
                          large.trees.p.acre <- ifelse(large.trees.p.acre >= 10, 1, 0)
                          
                          max.age <- max(c(raw.stand.age,local.trees$tree.age), na.rm=TRUE) #determine the age of the oldest tree in the stand
                          max.age <- ifelse(max.age >= 150, 1, 0)
                          
                          local.MOG.status <- ifelse(sum(large.trees.p.acre, max.age, na.rm=TRUE) == 2, 1, 0) #assign old growth status (1) if threshold is met
                          MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                          
                        }#end sub-sub-assessment for if site class is low or medium
                      }#end sub-assessment for white fir/grand fir in central part of region
                      
                      if(regional.geog %in% "non-central"){#begin sub-assessment for white fir/grand fir in non-central part of region
                        
                        if(site.class %in% c(1,2)){#begin sub-sub-assessment for if site class is high
                          
                          local.large.trees <- local.trees[local.trees$DIA >= 21, ] #isolate to large trees specific to regional and forest type threshold
                          large.trees.p.acre <- nrow(local.large.trees)/condition.area #convert from count to n per acre
                          large.trees.p.acre <- ifelse(large.trees.p.acre >= 20, 1, 0)
                          
                          max.age <- max(c(raw.stand.age,local.trees$tree.age), na.rm=TRUE) #determine the age of the oldest tree in the stand
                          max.age <- ifelse(max.age >= 150, 1, 0)
                          
                          local.MOG.status <- ifelse(sum(large.trees.p.acre, max.age, na.rm=TRUE) == 2, 1, 0) #assign old growth status (1) if threshold is met
                          MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                          
                        }#end sub-sub-assessment for if site class is high
                        
                        if(!(site.class %in% c(1,2))){#begin sub-sub-assessment for if site class is low or medium
                          
                          local.large.trees <- local.trees[local.trees$DIA >= 21, ] #isolate to large trees specific to regional and forest type threshold
                          large.trees.p.acre <- nrow(local.large.trees)/condition.area #convert from count to n per acre
                          large.trees.p.acre <- ifelse(large.trees.p.acre >= 10, 1, 0)
                          
                          max.age <- max(c(raw.stand.age,local.trees$tree.age), na.rm=TRUE) #determine the age of the oldest tree in the stand
                          max.age <- ifelse(max.age >= 150, 1, 0)
                          
                          local.MOG.status <- ifelse(sum(large.trees.p.acre, max.age, na.rm=TRUE) == 2, 1, 0) #assign old growth status (1) if threshold is met
                          MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                          
                        }#end sub-sub-assessment for if site class is low or medium
                      }#end sub-assessment for white fir/grand fir in non-central part of region
                      
                    }#end old growth instructions for (pacific northwest) white fir and grand fir forests outside of the management zone
                    
                    if(paz.group %in% "douglas fir"){#begin old growth instructions for (pacific northwest) douglas fir forests outside of the management zone
                      
                      #there is no site classification threshold for this forest type
                      
                      local.large.trees <- local.trees[local.trees$DIA >= 21, ] #isolate to large trees specific to regional and forest type threshold
                      large.trees.p.acre <- nrow(local.large.trees)/condition.area #convert from count to n per acre
                      large.trees.p.acre <- ifelse(large.trees.p.acre >= 8, 1, 0)
                      
                      max.age <- max(c(raw.stand.age,local.trees$tree.age), na.rm=TRUE) #determine the age of the oldest tree in the stand
                      max.age <- ifelse(max.age >= 150, 1, 0)
                      
                      local.MOG.status <- ifelse(sum(large.trees.p.acre, max.age, na.rm=TRUE) == 2, 1, 0) #assign old growth status (1) if threshold is met
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end old growth instructions for (pacific northwest) douglas fir forests outside of the management zone
                    
                    if(paz.group %in% "lodgepole pine"){#begin old growth instructions for (pacific northwest) lodgepole pine forests outside of the management zone
                      
                      #there is no site classification threshold for this forest type
                      
                      local.large.trees <- local.trees[local.trees$DIA >= 12, ] #isolate to large trees specific to regional and forest type threshold
                      large.trees.p.acre <- nrow(local.large.trees)/condition.area #convert from count to n per acre
                      large.trees.p.acre <- ifelse(large.trees.p.acre >= 60, 1, 0)
                      
                      max.age <- max(c(raw.stand.age,local.trees$tree.age), na.rm=TRUE) #determine the age of the oldest tree in the stand
                      max.age <- ifelse(max.age >= 120, 1, 0)
                      
                      local.MOG.status <- ifelse(sum(large.trees.p.acre, max.age, na.rm=TRUE) == 2, 1, 0) #assign old growth status (1) if threshold is met
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end old growth instructions for (pacific northwest) lodgepole pine forests outside of the management zone
                    
                    if(paz.group %in% "silver fir"){#begin old growth instructions for (pacific northwest) silver fir forests outside of the management zone
                      
                      if(site.class %in% 5){#begin sub-assessment for silver fir site class 5
                        
                        local.large.trees <- local.trees[local.trees$DIA >= 22, ] #isolate to large trees specific to regional and forest type threshold
                        large.trees.p.acre <- nrow(local.large.trees)/condition.area #convert from count to n per acre
                        large.trees.p.acre <- ifelse(large.trees.p.acre >= 9, 1, 0)
                        
                        max.age <- max(c(raw.stand.age,local.trees$tree.age), na.rm=TRUE) #determine the age of the oldest tree in the stand
                        max.age <- ifelse(max.age >= 260, 1, 0)
                        
                        local.MOG.status <- ifelse(sum(large.trees.p.acre, max.age, na.rm=TRUE) == 2, 1, 0) #assign old growth status (1) if threshold is met
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      }#end sub-assessment for silver fir site class 5
                      
                      if(site.class %in% 6){#begin sub-assessment for silver fir site class 6
                        
                        local.large.trees <- local.trees[local.trees$DIA >= 22, ] #isolate to large trees specific to regional and forest type threshold
                        large.trees.p.acre <- nrow(local.large.trees)/condition.area #convert from count to n per acre
                        large.trees.p.acre <- ifelse(large.trees.p.acre >= 1, 1, 0)
                        
                        max.age <- max(c(raw.stand.age,local.trees$tree.age), na.rm=TRUE) #determine the age of the oldest tree in the stand
                        max.age <- ifelse(max.age >= 360, 1, 0)
                        
                        local.MOG.status <- ifelse(sum(large.trees.p.acre, max.age, na.rm=TRUE) == 2, 1, 0) #assign old growth status (1) if threshold is met
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      }#end sub-assessment for silver fir site class 6
                      
                      if(site.class %in% c(2,3)){#begin sub-assessment for silver fir site classes 2 and 3
                        
                        local.large.trees <- local.trees[local.trees$DIA >= 26, ] #isolate to large trees specific to regional and forest type threshold
                        large.trees.p.acre <- nrow(local.large.trees)/condition.area #convert from count to n per acre
                        large.trees.p.acre <- ifelse(large.trees.p.acre >= 6, 1, 0)
                        
                        max.age <- max(c(raw.stand.age,local.trees$tree.age), na.rm=TRUE) #determine the age of the oldest tree in the stand
                        max.age <- ifelse(max.age >= 180, 1, 0)
                        
                        local.MOG.status <- ifelse(sum(large.trees.p.acre, max.age, na.rm=TRUE) == 2, 1, 0) #assign old growth status (1) if threshold is met
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      }#end sub-assessment for silver fir site classes 2 and 3
                      
                      if(site.class %in% 4){#begin sub-assessment for silver fir site class 4
                        
                        local.large.trees <- local.trees[local.trees$DIA >= 25, ] #isolate to large trees specific to regional and forest type threshold
                        large.trees.p.acre <- nrow(local.large.trees)/condition.area #convert from count to n per acre
                        large.trees.p.acre <- ifelse(large.trees.p.acre >= 7, 1, 0)
                        
                        max.age <- max(c(raw.stand.age,local.trees$tree.age), na.rm=TRUE) #determine the age of the oldest tree in the stand
                        max.age <- ifelse(max.age >= 200, 1, 0)
                        
                        local.MOG.status <- ifelse(sum(large.trees.p.acre, max.age, na.rm=TRUE) == 2, 1, 0) #assign old growth status (1) if threshold is met
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      }#end sub-assessment for silver fir site class 4
                    }#end old growth instructions for (pacific northwest) silver fir forests outside of the management zone
                    
                    
                    if(paz.group %in% "ponderosa pine"){#begin old growth instructions for (pacific northwest) ponderosa pine forests outside of the management zone
                      
                      if(site.class %in% c(1,2,3,4)){#begin sub-assessment for ponderosa pine classes 1,2,3,4 (med and high)
                        
                        #note that it is unclear what is meany by "very late decadent" in Table 14. Thus, both sets are thresholds are included here
                        
                        #normal
                        local.large.trees <- local.trees[local.trees$DIA >= 21, ] #isolate to large trees specific to regional and forest type threshold
                        large.trees.p.acre <- nrow(local.large.trees)/condition.area #convert from count to n per acre
                        large.trees.p.acre <- ifelse(large.trees.p.acre >= 13, 1, 0)
                        
                        max.age <- max(c(raw.stand.age,local.trees$tree.age), na.rm=TRUE) #determine the age of the oldest tree in the stand
                        max.age <- ifelse(max.age >= 150, 1, 0)
                        
                        local.MOG.status <- ifelse(sum(large.trees.p.acre, max.age, na.rm=TRUE) == 2, 1, 0) #assign old growth status (1) if threshold is met
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                        #very late decadent
                        local.large.trees <- local.trees[local.trees$DIA >= 31, ] #isolate to large trees specific to regional and forest type threshold
                        large.trees.p.acre <- nrow(local.large.trees)/condition.area #convert from count to n per acre
                        large.trees.p.acre <- ifelse(large.trees.p.acre >= 3, 1, 0)
                        
                        max.age <- max(c(raw.stand.age,local.trees$tree.age), na.rm=TRUE) #determine the age of the oldest tree in the stand
                        max.age <- ifelse(max.age >= 200, 1, 0)
                        
                        local.MOG.status <- ifelse(sum(large.trees.p.acre, max.age, na.rm=TRUE) == 2, 1, 0) #assign old growth status (1) if threshold is met
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end sub-assessment for ponderosa pine classes 1,2,3,4 (med and high)
                      
                      if(site.class > 4){#begin sub-assessment for ponderosa pine classes 4+ (low)
                        
                        #note that it is unclear what is meany by "very late decadent" in Table 14. Thus, both sets are thresholds are included here
                        
                        #normal
                        local.large.trees <- local.trees[local.trees$DIA >= 21, ] #isolate to large trees specific to regional and forest type threshold
                        large.trees.p.acre <- nrow(local.large.trees)/condition.area #convert from count to n per acre
                        large.trees.p.acre <- ifelse(large.trees.p.acre >= 10, 1, 0)
                        
                        max.age <- max(c(raw.stand.age,local.trees$tree.age), na.rm=TRUE) #determine the age of the oldest tree in the stand
                        max.age <- ifelse(max.age >= 150, 1, 0)
                        
                        local.MOG.status <- ifelse(sum(large.trees.p.acre, max.age, na.rm=TRUE) == 2, 1, 0) #assign old growth status (1) if threshold is met
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                        #very late decadent
                        local.large.trees <- local.trees[local.trees$DIA >= 31, ] #isolate to large trees specific to regional and forest type threshold
                        large.trees.p.acre <- nrow(local.large.trees)/condition.area #convert from count to n per acre
                        large.trees.p.acre <- ifelse(large.trees.p.acre >= 2, 1, 0)
                        
                        max.age <- max(c(raw.stand.age,local.trees$tree.age), na.rm=TRUE) #determine the age of the oldest tree in the stand
                        max.age <- ifelse(max.age >= 200, 1, 0)
                        
                        local.MOG.status <- ifelse(sum(large.trees.p.acre, max.age, na.rm=TRUE) == 2, 1, 0) #assign old growth status (1) if threshold is met
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end sub-assessment for ponderosa pine classes 4+ (low)
                    }#end old growth instructions for (pacific northwest) ponderosa pine forests outside of the management zone
                    
                    if(paz.group %in% "subalpine fir"){#begin old growth instructions for (pacific northwest) subalpine fir forests outside of the management zone
                      
                      #note that there are no instructions for subalpine fir with medium site class
                      
                      if(site.class %in% c(1,2)){#begin sub-assessment for subalpine fir in classes 1 and 2 (high)
                        
                        local.large.trees <- local.trees[local.trees$DIA >= 21, ] #isolate to large trees specific to regional and forest type threshold
                        large.trees.p.acre <- nrow(local.large.trees)/condition.area #convert from count to n per acre
                        large.trees.p.acre <- ifelse(large.trees.p.acre >= 10, 1, 0)
                        
                        max.age <- max(c(raw.stand.age,local.trees$tree.age), na.rm=TRUE) #determine the age of the oldest tree in the stand
                        max.age <- ifelse(max.age >= 150, 1, 0)
                        
                        local.MOG.status <- ifelse(sum(large.trees.p.acre, max.age, na.rm=TRUE) == 2, 1, 0) #assign old growth status (1) if threshold is met
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end sub-assessment for subalpine fir in classes 1 and 2 (high)
                      
                      if(site.class > 4){#begin sub-assessment for subalpine fir in classes 4+ (low)
                        
                        local.large.trees <- local.trees[local.trees$DIA >= 13, ] #isolate to large trees specific to regional and forest type threshold
                        large.trees.p.acre <- nrow(local.large.trees)/condition.area #convert from count to n per acre
                        large.trees.p.acre <- ifelse(large.trees.p.acre >= 10, 1, 0)
                        
                        max.age <- max(c(raw.stand.age,local.trees$tree.age), na.rm=TRUE) #determine the age of the oldest tree in the stand
                        max.age <- ifelse(max.age >= 150, 1, 0)
                        
                        local.MOG.status <- ifelse(sum(large.trees.p.acre, max.age, na.rm=TRUE) == 2, 1, 0) #assign old growth status (1) if threshold is met
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end sub-assessment for subalpine fir in classes 4+ (low)
                    }#end old growth instructions for (pacific northwest) subalpine fir forests outside of the management zone
                    
                    if(paz.group %in% "western hemlock"){#begin old growth instructions for (pacific northwest) western hemlock forests outside of the management zone
                      
                      if(site.class %in% 1){#begin sub-assessment for western hemlock in site class 1
                        
                        local.large.trees <- local.trees[local.trees$DIA >= 42, ] #isolate to large trees specific to regional and forest type threshold
                        large.trees.p.acre <- nrow(local.large.trees)/condition.area #convert from count to n per acre
                        large.trees.p.acre <- ifelse(large.trees.p.acre >= 8, 1, 0)
                        
                        max.age <- max(c(raw.stand.age,local.trees$tree.age), na.rm=TRUE) #determine the age of the oldest tree in the stand
                        max.age <- ifelse(max.age >= 200, 1, 0)
                        
                        local.MOG.status <- ifelse(sum(large.trees.p.acre, max.age, na.rm=TRUE) == 2, 1, 0) #assign old growth status (1) if threshold is met
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end sub-assessment for western hemlock in site class 1
                      
                      if(site.class %in% 2){#begin sub-assessment for western hemlock in site class 2
                        
                        local.large.trees <- local.trees[local.trees$DIA >= 35, ] #isolate to large trees specific to regional and forest type threshold
                        large.trees.p.acre <- nrow(local.large.trees)/condition.area #convert from count to n per acre
                        large.trees.p.acre <- ifelse(large.trees.p.acre >= 8, 1, 0)
                        
                        max.age <- max(c(raw.stand.age,local.trees$tree.age), na.rm=TRUE) #determine the age of the oldest tree in the stand
                        max.age <- ifelse(max.age >= 200, 1, 0)
                        
                        local.MOG.status <- ifelse(sum(large.trees.p.acre, max.age, na.rm=TRUE) == 2, 1, 0) #assign old growth status (1) if threshold is met
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end sub-assessment for western hemlock in site class 2
                      
                      if(site.class %in% 3){#begin sub-assessment for western hemlock in site class 3
                        
                        local.large.trees <- local.trees[local.trees$DIA >= 31, ] #isolate to large trees specific to regional and forest type threshold
                        large.trees.p.acre <- nrow(local.large.trees)/condition.area #convert from count to n per acre
                        large.trees.p.acre <- ifelse(large.trees.p.acre >= 8, 1, 0)
                        
                        max.age <- max(c(raw.stand.age,local.trees$tree.age), na.rm=TRUE) #determine the age of the oldest tree in the stand
                        max.age <- ifelse(max.age >= 200, 1, 0)
                        
                        local.MOG.status <- ifelse(sum(large.trees.p.acre, max.age, na.rm=TRUE) == 2, 1, 0) #assign old growth status (1) if threshold is met
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end sub-assessment for western hemlock in site class 3
                      
                      if(site.class %in% c(4,5)){#begin sub-assessment for western hemlock in site classes 4 and 5
                        
                        local.large.trees <- local.trees[local.trees$DIA >= 21, ] #isolate to large trees specific to regional and forest type threshold
                        large.trees.p.acre <- nrow(local.large.trees)/condition.area #convert from count to n per acre
                        large.trees.p.acre <- ifelse(large.trees.p.acre >= 8, 1, 0)
                        
                        max.age <- max(c(raw.stand.age,local.trees$tree.age), na.rm=TRUE) #determine the age of the oldest tree in the stand
                        max.age <- ifelse(max.age >= 200, 1, 0)
                        
                        local.MOG.status <- ifelse(sum(large.trees.p.acre, max.age, na.rm=TRUE) == 2, 1, 0) #assign old growth status (1) if threshold is met
                        MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                        
                      }#end sub-assessment for western hemlock in site classes 4 and 5
                    }#end old growth instructions for (pacific northwest) western hemlock forests outside of the management zone
                  }#end instructions for if the plot is outside the northwest forest plan
                  
                  #if the plot is not old growth, test to see if it is at least mature based on Table 19 of MOG document
                  if(local.MOG.status == 0){#begin (pacific northwest) mature assessment
                    
                    if(forest.type %in% c(700:722,900:905,960:976,920:935)){#begin mature assessment for (pacific northwest) hardwoods
                      
                      t.ddiscore <- ifelse(ddiscore >= 47.5, 0.31, 0)
                      t.HTquart <- ifelse(HTquart >= 68.5, 0.31, 0)
                      t.badom <- ifelse(badom >= 150, 0.21, 0)
                      t.snagbatot <- ifelse(snagbatot >= 24.9, 0.16, 0)
                      
                      local.MOG.status <- sum(t.ddiscore,t.HTquart,t.badom,t.snagbatot,na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for (pacific northwest) hardwoods
                    
                    if(forest.type %in% c(360:369, 180:185)){#begin mature assessment for (pacific northwest) softwoods
                      
                      t.QMDdom <- ifelse(QMDdom >= 14.2, 0.54, 0)
                      t.badom <- ifelse(badom >= 30.9, 0.46, 0)
                      
                      local.MOG.status <- sum(t.QMDdom,t.badom,na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for (pacific northwest) softwoods
                    
                    if(paz.group %in% c("douglas fir") & nrow(nwfp.test) == 0){#begin mature assessment for (pacific northwest) interior douglas fir
                      
                      t.QMDdom <- ifelse(QMDdom >= 11.1, 0.42, 0)
                      t.ddiscore <- ifelse(ddiscore >= 30.2, 0.38, 0)
                      t.badom <- ifelse(badom >= 60.1, 0.21, 0)
                      
                      local.MOG.status <- sum(t.QMDdom,t.ddiscore,t.badom,na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for (pacific northwest) interior douglas fir
                    
                    if(forest.type %in% c(200:203)){#begin mature assessment for (pacific northwest) douglas fir (regardless of locale)
                      
                      t.QMDdom <- ifelse(QMDdom >= 11.1, 0.42, 0)
                      t.ddiscore <- ifelse(ddiscore >= 30.2, 0.38, 0)
                      t.badom <- ifelse(badom >= 60.1, 0.21, 0)
                      
                      local.MOG.status <- sum(t.QMDdom,t.ddiscore,t.badom,na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for (pacific northwest) douglas fir (regardless of locale)
                    
                    if(paz.group %in% c("douglas fir") & nrow(nwfp.test) > 0){#begin mature assessment for (pacific northwest) pacific douglas fir
                      
                      t.QMDdom <- ifelse(QMDdom >= 12.7, 0.45, 0)
                      t.ddiscore <- ifelse(ddiscore >= 32.6, 0.33, 0)
                      t.badom <- ifelse(badom >= 42.3, 0.23, 0)
                      
                      local.MOG.status <- sum(t.QMDdom,t.ddiscore,t.badom,na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for (pacific northwest) pacific douglas fir
                    
                    if(paz.group %in% c("mountain hemlock") | forest.type %in% c(260:271)){#begin mature assessment for (pacific northwest) mountain hemlock
                      
                      t.QMDdom <- ifelse(QMDdom >= 13.1, 0.29, 0)
                      t.badom <- ifelse(badom >= 126.6, 0.2, 0)
                      t.HTsd <- ifelse(HTsd >= 58.5, 0.2, 0)
                      t.HTquart <- ifelse(HTquart >= 42.7, 0.19, 0)
                      t.tpadom <- ifelse(tpadom <= 77.4, 0.12, 0)
                      
                      local.MOG.status <- sum(t.QMDdom,t.badom,t.HTsd,t.HTquart,t.tpadom,na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for (pacific northwest) mountain hemlock
                    
                    if(paz.group %in% c("ponderosa pine","jeffrey pine - knobcone pine","lodgepole pine") | forest.type %in% c(220:226)){#begin mature assessment for (pacific northwest) ponderosa and lodgepole pine
                      
                      t.QMDdom <- ifelse(QMDdom >= 7.7, 0.34, 0)
                      t.ddiscore <- ifelse(ddiscore >= 15.3, 0.28, 0)
                      t.HTsd <- ifelse(HTsd >= 31.2, 0.22, 0)
                      t.tpadom <- ifelse(tpadom <= 31.7, 0.16, 0)
                      
                      local.MOG.status <- sum(t.QMDdom,t.ddiscore,t.HTsd,t.tpadom,na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for (pacific northwest) ponderosa and lodgepole pine
                    
                    #note that there is specific mature criteria for "very late decadent" ponderosa pines. No guidance is given as to what that means. Any plot that would meet that maturity critera would also meet the standard ponderosa critera. Thus, it is omitted here.
                    
                    if(paz.group %in% c("port orford cedar","redwood")){#begin mature assessment for (pacific northwest) port orford cedar and redwood
                      
                      t.ddiscore <- ifelse(ddiscore >= 44.4, 0.62, 0)
                      t.QMDdom <- ifelse(QMDdom >= 13, 0.38, 0)
                      
                      local.MOG.status <- sum(t.QMDdom,t.ddiscore,na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for (pacific northwest) port orford cedar and redwood
                    
                    if(paz.group %in% c("shasta red fir","silver fir")){#begin mature assessment for (pacific northwest) red and silver fir
                      
                      t.QMDdom <- ifelse(QMDdom >= 17.1, 0.29, 0)
                      t.HTsd <- ifelse(t.HTsd >= 72.2, 0.2, 0)
                      t.badom <- ifelse(badom >= 161.6, 0.19, 0)
                      t.snagbatot <- ifelse(snagbatot >= 39.7, 0.18, 0)
                      t.tpadom <- ifelse(tpadom <= 53.1, 0.14, 0)
                      
                      local.MOG.status <- sum(t.QMDdom,t.HTsd,t.badom,t.snagbatot,t.tpadom,na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for (pacific northwest) red and silver fir
                    
                    if(paz.group %in% c("sitka spruce")){#begin mature assessment for (pacific northwest) sitka spruce
                      
                      t.QMDdom <- ifelse(QMDdom >= 24.3, 0.3, 0)
                      t.HTsd <- ifelse(t.HTsd >= 63.5, 0.22, 0)
                      t.badom <- ifelse(badom >= 184.6, 0.2, 0)
                      t.snagbatot <- ifelse(snagbatot >= 54.5, 0.13, 0)
                      t.tpadom <- ifelse(tpadom <= 37.6, 0.15, 0)
                      
                      local.MOG.status <- sum(t.QMDdom,t.HTsd,t.badom,t.snagbatot,t.tpadom,na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for (pacific northwest) sitka spruce
                    
                    if(paz.group %in% c("subalpine fir")){#begin mature assessment for (pacific northwest) subalpine fir
                      
                      t.ddiscore <- ifelse(ddiscore >= 27.8, 0.4, 0)
                      t.HTquart <- ifelse(HTquart >= 39.8, 0.32, 0)
                      t.HTsd <- ifelse(HTsd >= 41.3, 0.29, 0)
                      
                      local.MOG.status <- sum(t.ddiscore,t.HTquart,t.HTsd,na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for (pacific northwest) subalpine fir
                    
                    
                    if(forest.type %in% c(265,266)){#begin mature assessment for (pacific northwest) engelmann spruce
                      
                      t.ddiscore <- ifelse(ddiscore >= 33.2, 0.42, 0)
                      t.QMDdom <- ifelse(QMDdom >= 8.8, 0.35, 0)
                      t.HTsd <- ifelse(HTsd >= 42.9, 0.23, 0)
                      
                      local.MOG.status <- sum(t.ddiscore,t.QMDdom,t.HTsd,na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for (pacific northwest) engelmann spruce
                    
                    if(paz.group %in% c("tanoak")){#begin mature assessment for (pacific northwest) tanoak
                      
                      t.QMDdom <- ifelse(QMDdom >= 15.3, 0.29, 0)
                      t.ddiscore <- ifelse(ddiscore >= 56, 0.24, 0)
                      t.HTquart <- ifelse(HTquart >= 51.7, 0.16, 0)
                      t.tpadom <- ifelse(tpadom <= 55.9, 0.16, 0)
                      t.HTsd <- ifelse(HTsd >= 64, 0.15, 0)
                      
                      local.MOG.status <- sum(t.QMDdom,t.ddiscore,t.HTquart,t.tpadom,t.HTsd,na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for (pacific northwest) tanoak
                    
                    if(paz.group %in% c("western hemlock")){#begin mature assessment for (pacific northwest) western hemlock
                      
                      t.QMDdom <- ifelse(QMDdom >= 19.9, 0.33, 0)
                      t.badom <- ifelse(badom >= 156.2, 0.2, 0)
                      t.HTsd <- ifelse(HTsd >= 25.9, 0.17, 0)
                      t.snagbatot <- ifelse(snagbatot >= 63.2, 0.17, 0)
                      t.tpadom <- ifelse(tpadom <= 42, 0.14, 0)
                      
                      local.MOG.status <- sum(t.QMDdom,t.badom,t.HTsd,t.snagbatot,t.tpadom,na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for (pacific northwest) western hemlock
                    
                    if(paz.group %in% c("white fir - grand fir") | forest.type %in% c(267,261)){#begin mature assessment for (pacific northwest) white and grand fir
                      
                      t.QMDdom <- ifelse(QMDdom >= 12.3, 0.33, 0)
                      t.ddiscore <- ifelse(ddiscore >= 40.1, 0.31, 0)
                      t.HTsd <- ifelse(HTsd >= 46.8, 0.2, 0)
                      t.snagbatot <- ifelse(snagbatot >= 8.6, 0.16, 0)
                      
                      local.MOG.status <- sum(t.QMDdom,t.ddiscore,t.HTsd,t.snagbatot,na.rm=TRUE)
                      MOG.vector <- c(MOG.vector,local.MOG.status) #add local MOG status to vector
                      
                    }#end mature assessment for (pacific northwest) white and grand fir
                    
                  }#end (pacific northwest) mature assessment
                } #end instructions for pacific northwest region
                
              }#end processing if all tree metrics are available
            }#end calculating indices/metrics if data are available
          } #end instructions for if there are trees in the condition
        }}#end looping through conditions (c) if conditions are present
      
      if(length(MOG.vector) > 1){FIA.plots[p,"p.mog"] <- max(MOG.vector, na.rm=TRUE)} #assign the highest MOG rating for this plot as the final MOG status
      if(length(MOG.vector) == 1){FIA.plots[p,"p.mog"] <- NA} #assign the highest MOG rating for this plot as the final MOG status
      
      if(p <= 10){stop.time <- Sys.time()}
      if(p <= 10){time.vector <- c(time.vector, difftime(stop.time,start.time,units = "min"))}
      
      if(p == 1){message("FIA data processing has begun - progress will print periodically")} #print initial message indicating R is processing the FIA plots
      
      #warn how long this will take
      if(p == 10){ #begin time warning
        
        per.point <- mean(as.numeric(time.vector),na.rm=TRUE) 
        est.time <- as.numeric(per.point*nrow(FIA.plots)) #multiply time difference by n points in the data set
        est.time <- est.time*1.25 #inflate estimated time since R consistently under estimates how long this will take
        
        message(paste("Approximately ", round(est.time), " minutes to completion (based on processing time of the first 10 FIA plots, actual processing time may vary).", sep=""))
        
      } #end time warning
      
      #print progress message
      if(p == floor(.1*nrow(FIA.plots))){message("10% of FIA plots have been processed")}
      if(p == floor(.2*nrow(FIA.plots))){message("20% of FIA plots have been processed")}
      if(p == floor(.3*nrow(FIA.plots))){message("30% of FIA plots have been processed")}
      if(p == floor(.4*nrow(FIA.plots))){message("40% of FIA plots have been processed")}
      if(p == floor(.5*nrow(FIA.plots))){message("50% of FIA plots have been processed")}
      if(p == floor(.6*nrow(FIA.plots))){message("60% of FIA plots have been processed")}
      if(p == floor(.7*nrow(FIA.plots))){message("70% of FIA plots have been processed")}
      if(p == floor(.8*nrow(FIA.plots))){message("80% of FIA plots have been processed")}
      if(p == floor(.9*nrow(FIA.plots))){message("90% of FIA plots have been processed")}
      if(p == nrow(FIA.plots)){message("100% of FIA plots have been processed")}
      
    }#end looping through plots (p)
  } #end instructions for if rFIA package is being used
  
  FIA.plots <- FIA.plots[!(is.na(FIA.plots$p.mog)),] #remove plots that have no MOG data (missing survey data)
  FIA.plots <- FIA.plots[FIA.plots$p.mog >= 0, ] #removes any plots that have negative MOG data (missing survey data)
  
  #leave for the user
  MOG.points <- FIA.plots[,c("p.mog","geometry")]
  names(MOG.points) <- c("p.mog","geometry")
  MOG.points <<- MOG.points
  
  if(interpolate == TRUE){ #begin interpolation process if specified by used (default = TRUE)
    
    if(max(MOG.points$p.mog)<0.5){stop("No mature or old growth forest detected. Imputation cancelled.")} #print error message if interpolation isn't worth it
    
    #print progress message
    message("Accessing landcover data for interpolation (interpolate = TRUE)")
    
    NLCD.small <- terra::rast(canopy.path) #load canopy data
    temp.locale.vect <- vect(st_transform(x= locale, #convert locale to match crs of NLCD and work with terra
                                          crs = st_crs(NLCD.small)))
    
    NLCD.small <- terra::crop(x = NLCD.small, #crop raster to locale
                              y = temp.locale.vect)
    
    NLCD.large <- terra::aggregate(x = NLCD.small,
                                   fact = (resolution/30), #specify the number of cells to aggregate horizontally and vertically. Default if 53 cells for a 1590 meter (~1 mile) grid (accommodating the 0.5 mile displacement of FIA plot data)
                                   fun = "mean", #takes the mean value of all the 30x30 m cells making up the new grid cell
                                   na.rm=TRUE) #removing NA's
    
    message("Calculating ecosystem and region-specific canopy cover threshold (interpolate = TRUE)")
    #assess regional canopy cover threshold
    forested.plots <- merge(x = FIA.plots, #merge conditions to plots
                            y = FIA.condition,
                            by = "PLOT")
    forested.plots <- forested.plots[forested.plots$COND_STATUS_CD == 1, ] #retain only conditions which are classified as forested land
    forested.plots <- forested.plots[match(unique(forested.plots$PLOT), forested.plots$PLOT), ] #retain only the first instance of each plot (as to avoid duplicates from plots with multiple conditions or surveyed over multiple years)
    
    forested.plots <- terra::vect(forested.plots) #convert plots to terra format
    forested.plots <- terra::project(x = forested.plots, #reproject plots to crs of NLCD raster
                                     y = crs(NLCD.small))
    
    ecoregions <- st_read(paste(source.path,"/utility_ecoregions/USA_ecoregions.shp", sep = ""))
    ecoregions <- st_transform(x = ecoregions,
                               crs = st_crs(locale)) #transform ecoregions to match existing data
    ecoregions <- st_filter(x = ecoregions, #isolate only ecoregions that overlap with forested points
                            y = locale)
    eco.types <- as.vector(unique(ecoregions$NA_L2NAME)) #create a vector of ecotypes in the locale
    eco.types <- eco.types[which(!(eco.types %in% "WATER"))] #remove water from potential ecotypes
    
    temp.locale <- terra::project(x = vect(locale), #reproject locale (faster than reporjecting raster)
                                  y = crs(NLCD.large))
    
    NLCD.trees.vect <- data.frame() #create an empty data frame to fill
    
    for(e in 1:length(eco.types)){#begin looping through ecoregions for region-specific interpolation
      
      local.ecoregion <- ecoregions[ecoregions$NA_L2NAME %in% eco.types[e], ] #isolate to local ecoregion type
      local.ecoregion <- st_union(local.ecoregion) #dissolve to one shapefile
      local.ecoregion <- st_as_sf(local.ecoregion) #convert back to shapefile format
      local.ecoregion <- st_transform(x = local.ecoregion, #transform ecoregion to match forested plot transormation
                                      crs = st_crs(forested.plots))
      local.forested.plots <- st_filter(x = st_as_sf(forested.plots), #isolate to only forested points in the specific ecoregion
                                        y = local.ecoregion)
      
      canopy.threshold <- NA #assign NA now to resolve downstream issues
      if(nrow(local.forested.plots)>0){
        
        canopy.threshold <- terra::extract(x = NLCD.small, #extract canopy value at each unique forested FIA plot within local ecoregion
                                           y = vect(local.forested.plots))
        canopy.threshold <- mean(as.numeric(canopy.threshold[,2]), na.rm=TRUE) #find mean percent canopy cover to use as a threshold
        canopy.threshold <- ifelse(canopy.threshold >= 50, 50, canopy.threshold)} #if mean percent canopy is greater than 50, drop the threshold to 50% to capture more forest (the larger concern is finding an adequate threshold when forest is scarce, not when it is abundant)
      
      if(is.na(canopy.threshold)){canopy.threshold <- 200} #if there is absolutely no canopy cover, all values are NA so mean is NA. Set cutoff to 200% to not allow any data in this region to be interpolated as forest.
      
      NLCD <- mask(x = NLCD.large,
                   mask = local.ecoregion)
      
      reclass.rules <- as.data.frame(matrix(c(canopy.threshold,100))) #save only cells with at least 50% forest cover but no more than 100% (values over water recieve erroneous values greater than 100)
      
      #reclassify raster based on previously established rules
      NLCD.trees <- terra::classify(x = NLCD, #raster to reclassify
                                    rcl = reclass.rules) #rules for how to reclassify it
      
      #convert raster to points
      local.NLCD.trees.vect <- terra::as.points(NLCD.trees) #convert tree raster to points
      local.NLCD.trees.vect <- st_as_sf(local.NLCD.trees.vect) #convert back to sf for ease of saving
      if(nrow(local.NLCD.trees.vect)>0){NLCD.trees.vect <- rbind(NLCD.trees.vect,local.NLCD.trees.vect)} #add local ecoregion trees to locale trees (if data exist)
      
    } #end looping through ecoregions for region-specific interpolation
    
    message("Preparing landcover data for interpolation (interpolate = TRUE)")
    
    NLCD.trees.vect <- vect(NLCD.trees.vect) #convert back to terra
    
    NLCD.trees.vect <- terra::mask(x = NLCD.trees.vect, #clip tree raster points to the locale
                                   mask = temp.locale)
    NLCD.trees.vect$imputed <- NA
    
    FIA.plots.vect <- terra::vect(FIA.plots) #convert FIA plots to terra format
    FIA.plots.vect <- terra::project(x = FIA.plots.vect, #reproject FIA plots to tree points' crs
                                     y = crs(NLCD.trees.vect))
    
    nearby.plots <- as.data.frame(terra::nearby(x = NLCD.trees.vect, #find the nearest three FIA plots to each point
                                                y = FIA.plots.vect,
                                                k = 3))
    
    message("Predicting MOG status across locale (interpolate = TRUE)")
    
    #process each cell
    for(i in 1:nrow(nearby.plots)){ #begin looping through tree cell points
      
      if(i == 1){time.vector <- vector()}
      if(i <= 10){start.time <- Sys.time()}
      
      #calculate distance between tree raster point and nearest three FIA plots
      dist.k1 <- terra::distance(x = NLCD.trees.vect[i,],
                                 y = FIA.plots.vect[nearby.plots[i,"k1"]])[1,1]
      dist.k2 <- terra::distance(x = NLCD.trees.vect[i,],
                                 y = FIA.plots.vect[nearby.plots[i,"k2"]])[1,1]
      dist.k3 <- terra::distance(x = NLCD.trees.vect[i,],
                                 y = FIA.plots.vect[nearby.plots[i,"k3"]])[1,1]
      
      #determine MOG status of nearest three FIA plots
      status.k1 <- st_drop_geometry(st_as_sf(FIA.plots.vect[nearby.plots[i,"k1"],"p.mog"]))[1,1]
      status.k2 <- st_drop_geometry(st_as_sf(FIA.plots.vect[nearby.plots[i,"k2"],"p.mog"]))[1,1]
      status.k3 <- st_drop_geometry(st_as_sf(FIA.plots.vect[nearby.plots[i,"k3"],"p.mog"]))[1,1]
      
      #assign MOG value if cell contains an FIA plot
      if(dist.k1 <= res(NLCD.trees)[1]){NLCD.trees.vect[i, "imputed"] <- status.k1}
      if(dist.k2 <= res(NLCD.trees)[1]){NLCD.trees.vect[i, "imputed"] <- status.k2}
      if(dist.k3 <= res(NLCD.trees)[1]){NLCD.trees.vect[i, "imputed"] <- status.k3}
      if(dist.k1 <= res(NLCD.trees)[1] & dist.k2 <= res(NLCD.trees)[1]){NLCD.trees.vect[i, "imputed"] <- max(c(status.k1,status.k2))}
      if(dist.k1 <= res(NLCD.trees)[1] & dist.k3 <= res(NLCD.trees)[1]){NLCD.trees.vect[i, "imputed"] <- max(c(status.k1,status.k3))}
      if(dist.k2 <= res(NLCD.trees)[1] & dist.k3 <= res(NLCD.trees)[1]){NLCD.trees.vect[i, "imputed"] <- max(c(status.k2,status.k3))}
      if(dist.k1 <= res(NLCD.trees)[1] & dist.k2 <= res(NLCD.trees)[1] & dist.k3 <= res(NLCD.trees)[1]){NLCD.trees.vect[i, "imputed"] <- max(c(status.k1,status.k2,status.k3))}
      
      #assign MOG value if cell DOES NOT contain an FIA plot and needs imputation
      if(dist.k1 > res(NLCD.trees)[1] & dist.k2 > res(NLCD.trees)[1] & dist.k3 > res(NLCD.trees)[1]){ #begin conditions for if tree cell is not overlapping with an FIA plot
        #perform inverse distance weighting
        numerator <- (status.k1/ (dist.k1^power)) + (status.k2/ (dist.k2^power)) + (status.k3/ (dist.k1^power)) #specify power to make edges sharper or more gradual (higher values weight distance heavier)
        denomenator <- (1/(dist.k1^power)) + (1/(dist.k2^power)) + (1/(dist.k3^power)) #specify power to make edges sharper or more gradual (higher values weight distance heavier)
        NLCD.trees.vect[i, "imputed"] <- numerator/denomenator #assign imputed MOG value to cell
      } #end conditions for if tree cell is not overlapping with an FIA plot
      
      if(i <= 10){stop.time <- Sys.time()}
      if(i <= 10){time.vector <- c(time.vector, difftime(stop.time,start.time,units = "min"))}
      
      #warn how long this will take
      if(i == 10){ #begin time warning
        
        per.point <- mean(as.numeric(time.vector))
        est.time <- as.numeric(per.point*nrow(NLCD.trees.vect)) #multiply time difference by n points in the dataset
        est.time <- est.time*1.25 #inflate estimated time since R consistently underestimates time to completion
        
        message(paste("Interpolation has begun. Approximately ", round(est.time), " minutes to completion (based on processing time of the first 10 FIA plots, actual processing time may vary).", sep=""))
        
      } #end time warning
      
      if(i == floor(.1*nrow(NLCD.trees.vect))){message("10% of forested area has been interpolated")}
      if(i == floor(.2*nrow(NLCD.trees.vect))){message("20% of forested area has been interpolated")}
      if(i == floor(.3*nrow(NLCD.trees.vect))){message("30% of forested area has been interpolated")}
      if(i == floor(.4*nrow(NLCD.trees.vect))){message("40% of forested area has been interpolated")}
      if(i == floor(.5*nrow(NLCD.trees.vect))){message("50% of forested area has been interpolated")}
      if(i == floor(.6*nrow(NLCD.trees.vect))){message("60% of forested area has been interpolated")}
      if(i == floor(.7*nrow(NLCD.trees.vect))){message("70% of forested area has been interpolated")}
      if(i == floor(.8*nrow(NLCD.trees.vect))){message("80% of forested area has been interpolated")}
      if(i == floor(.9*nrow(NLCD.trees.vect))){message("90% of forested area has been interpolated")}
      if(i == nrow(NLCD.trees.vect)){message("100% of forested area has been interpolated")}
      
    } #end looping through tree cell points
    
    #Create raster
    raster.MOG <- terra::rasterize(x = NLCD.trees.vect, #input shapefile
                                   y = NLCD.trees, #raster to use as a template
                                   fun = "max",
                                   field = "imputed", #variable of x to use to populate cells
                                   touches = TRUE, #instructs R to consider any point that touches a particular cell
                                   background = 0) #assign any cell that doesn't have MOG data a value of zero (not mature/old growth)
    
    
    #clip raster to locale and save
    raster.MOG <- mask(x = raster.MOG,
                       mask = temp.locale)
    raster.MOG <- terra::project(x = raster.MOG, #SpatRaster to project
                                 y = crs(vect(locale))) #crs to project to
    names(raster.MOG) <- "MOG" #rename
    
    raster.MOG[raster.MOG > 1] <- 1 #if any cells go greater than 1 in the imputation process, bring them back down to 1 (old growth)
    
    #leave raster for user
    MOG.raster <<- raster.MOG
    
  } #end interpolation process if specified by used (default = TRUE)
  
  message("process complete - check environment for saved objects")
  if(length(northern.trigger)>0){message("please note that the classification of forests in the Northern Region using FIA data is an approximation - a field assessment is required for precise classification")}
  
} #end function

getOlympic <- function(source.path = NULL){ #begin function
  
  #throw error if source path is not provided
  if(is.null(source.path)){stop("source.path argument must be provided. Please provide the file path (in quotations) to the mapMog folder to access utility data required to execute the function.")}
  
  #load/install necessary packages
  if (!require(sf)) message("Installing required package 'sf.'")
  if (!require(sf)) install.packages('sf') #Does user have the package downloaded? If so move on, if not, download.
  library(sf) #load library
  
  olympicNF <- st_read(paste(source.path,"/tutorial_Olympic/OlympicNF.shp",sep="")) #read in sf data for Olympic National Park
  
  olympicNF #leave shapefile here for R to grab
  
} #end function
