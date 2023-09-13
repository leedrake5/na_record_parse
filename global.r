library(dplyr)
library(rgdal)
library(lipdR)

lm_eqn = function(m) {
    
    l <- list(a = format(coef(m)[1], digits = 2),
    b = format(abs(coef(m)[2]), digits = 2),
    r2 = format(summary(m)$r.squared, digits = 3));
    
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
    
    
    as.character(as.expression(eq));
}

my.max <- function(x) ifelse( !all(is.na(x)), max(x, na.rm=T), NA)
my.min <- function(x) ifelse( !all(is.na(x)), min(x, na.rm=T), NA)


###Dependent Transformation
scaleTransform <- function(values, the_min=NULL, the_max=NULL){
    
    if(is.null(the_min)){
        the_min <- my.min(values)
    }
    
    if(is.null(the_max)){
        the_max <- my.max(values)
    }
    
    train_scale <- ((values-the_min)/(the_max-the_min))

    return(train_scale)
}

scaleDecode <- function(values, y_min, y_max){

    (y_max-y_min)*values + y_min
    
    #y_train_decoded <- (values*(y_max-y_min)) + y_min

    #return(y_train_decoded)
}



load("data/wNAm1_0_0.RData")
metadata <- read.csv("data/wNAm_v1.csv")
metadata$ID <- formatC(make.names(metadata$siteName), width=5, format="d", flag="0")
archive_types <- unique(sapply(TS, function(x) x$archiveType))
proxy_types <- unique(unlist(sapply(TS, function(x) x$paleoData_proxyGeneral)))
site_names <- unique(lipdR::pullTsVariable(sTS,"geo_siteName"))
season_types <-  unique(lipdR::pullTsVariable(sTS,"interpretation1_seasonalityGeneral"))
variable_types <-  unique(lipdR::pullTsVariable(sTS,"paleoData_variableName"))

too_many_c <- c("deg C", "deg c", "Degrees Celcius", "degc", "Degrees celcius", "DegC")

varsColor <- c("archiveType", "proxy", "variableName", "seasonGeneral")

ug <- pullTsVariable(TS,"paleoData_units")
sg <- pullTsVariable(TS,"interpretation1_seasonalityGeneral")
pg <- pullTsVariable(TS,"paleoData_proxyGeneral")
ag <- pullTsVariable(TS,"archiveType")
doi <- pullTsVariable(TS,"pub1_doi")

unit_types <- unique(ug)
archive_types <- unique(ag)
season_types <- unique(sg)
proxy_types <- unique(pg)

mis_assign <- function(simple_merge){
    simple_merge$age <- simple_merge$age*-1000
    simple_merge$MIS_stage <- NA
    simple_merge[simple_merge$age > -14000 & simple_merge$age <= 0, "MIS_stage"] <- "MIS_001"
    simple_merge[simple_merge$age > -29000 & simple_merge$age <= -14000, "MIS_stage"] <- "MIS_002"
    simple_merge[simple_merge$age > -57000 & simple_merge$age <= -29000, "MIS_stage"] <- "MIS_003"
    simple_merge[simple_merge$age > -71000 & simple_merge$age <= -57000, "MIS_stage"] <- "MIS_004"
    simple_merge[simple_merge$age > -130000 & simple_merge$age <= -71000, "MIS_stage"] <- "MIS_005"
    simple_merge[simple_merge$age > -191000 & simple_merge$age <= -130000, "MIS_stage"] <- "MIS_006"
    simple_merge[simple_merge$age > -243000 & simple_merge$age <= -191000, "MIS_stage"] <- "MIS_007"
    simple_merge[simple_merge$age > -300000 & simple_merge$age <= -243000, "MIS_stage"] <- "MIS_008"
    simple_merge[simple_merge$age > -337000 & simple_merge$age <= -300000, "MIS_stage"] <- "MIS_009"
    simple_merge[simple_merge$age > -374000 & simple_merge$age <= -337000, "MIS_stage"] <- "MIS_010"
    simple_merge[simple_merge$age > -424000 & simple_merge$age <= -374000, "MIS_stage"] <- "MIS_011"
    simple_merge[simple_merge$age > -478000 & simple_merge$age <= -424000, "MIS_stage"] <- "MIS_012"
    simple_merge[simple_merge$age > -524000 & simple_merge$age <= -478000, "MIS_stage"] <- "MIS_013"
    simple_merge[simple_merge$age > -563000 & simple_merge$age <= -524000, "MIS_stage"] <- "MIS_014"
    simple_merge[simple_merge$age > -621000 & simple_merge$age <= -563000, "MIS_stage"] <- "MIS_015"
    simple_merge[simple_merge$age > -676000 & simple_merge$age <= -621000, "MIS_stage"] <- "MIS_016"
    simple_merge[simple_merge$age > -712000 & simple_merge$age <= -676000, "MIS_stage"] <- "MIS_017"
    simple_merge[simple_merge$age > -761000 & simple_merge$age <= -712000, "MIS_stage"] <- "MIS_018"
    simple_merge[simple_merge$age > -790000 & simple_merge$age <= -761000, "MIS_stage"] <- "MIS_019"
    simple_merge[simple_merge$age > -814000 & simple_merge$age <= -790000, "MIS_stage"] <- "MIS_020"
    simple_merge[simple_merge$age > -866000 & simple_merge$age <= -814000, "MIS_stage"] <- "MIS_021"
    simple_merge[simple_merge$age > -900000 & simple_merge$age <= -866000, "MIS_stage"] <- "MIS_022"
    simple_merge[simple_merge$age > -917000 & simple_merge$age <= -900000, "MIS_stage"] <- "MIS_023"
    simple_merge[simple_merge$age > -936000 & simple_merge$age <= -917000, "MIS_stage"] <- "MIS_024"
    simple_merge[simple_merge$age > -959000 & simple_merge$age <= -936000, "MIS_stage"] <- "MIS_025"
    simple_merge[simple_merge$age > -970000 & simple_merge$age <= -959000, "MIS_stage"] <- "MIS_026"
    simple_merge[simple_merge$age > -982000 & simple_merge$age <= -970000, "MIS_stage"] <- "MIS_027"
    simple_merge[simple_merge$age > -1014000 & simple_merge$age <= -982000, "MIS_stage"] <- "MIS_028"
    simple_merge[simple_merge$age > -1031000 & simple_merge$age <= -1014000, "MIS_stage"] <- "MIS_029"
    simple_merge[simple_merge$age > -1062000 & simple_merge$age <= -1031000, "MIS_stage"] <- "MIS_030"
    simple_merge[simple_merge$age > -1081000 & simple_merge$age <= -1062000, "MIS_stage"] <- "MIS_031"
    simple_merge[simple_merge$age > -1104000 & simple_merge$age <= -1081000, "MIS_stage"] <- "MIS_032"
    simple_merge[simple_merge$age > -1114000 & simple_merge$age <= -1104000, "MIS_stage"] <- "MIS_033"
    simple_merge[simple_merge$age > -1141000 & simple_merge$age <= -1114000, "MIS_stage"] <- "MIS_034"
    simple_merge[simple_merge$age > -1190000 & simple_merge$age <= -1141000, "MIS_stage"] <- "MIS_035"
    simple_merge[simple_merge$age > -1215000 & simple_merge$age <= -1190000, "MIS_stage"] <- "MIS_036"
    simple_merge[simple_merge$age > -1244000 & simple_merge$age <= -1215000, "MIS_stage"] <- "MIS_037"
    simple_merge[simple_merge$age > -1264000 & simple_merge$age <= -1244000, "MIS_stage"] <- "MIS_038"
    simple_merge[simple_merge$age > -1286000 & simple_merge$age <= -1264000, "MIS_stage"] <- "MIS_039"
    simple_merge[simple_merge$age > -1304000 & simple_merge$age <= -1286000, "MIS_stage"] <- "MIS_040"
    simple_merge[simple_merge$age > -1320000 & simple_merge$age <= -1304000, "MIS_stage"] <- "MIS_041"
    simple_merge[simple_merge$age > -1344000 & simple_merge$age <= -1320000, "MIS_stage"] <- "MIS_042"
    simple_merge[simple_merge$age > -1362000 & simple_merge$age <= -1344000, "MIS_stage"] <- "MIS_043"
    simple_merge[simple_merge$age > -1383000 & simple_merge$age <= -1362000, "MIS_stage"] <- "MIS_044"
    simple_merge[simple_merge$age > -1405000 & simple_merge$age <= -1383000, "MIS_stage"] <- "MIS_045"
    simple_merge[simple_merge$age > -1424000 & simple_merge$age <= -1405000, "MIS_stage"] <- "MIS_046"
    simple_merge[simple_merge$age > -1452000 & simple_merge$age <= -1424000, "MIS_stage"] <- "MIS_047"
    simple_merge[simple_merge$age > -1469000 & simple_merge$age <= -1452000, "MIS_stage"] <- "MIS_048"
    simple_merge[simple_merge$age > -1492000 & simple_merge$age <= -1469000, "MIS_stage"] <- "MIS_049"
    simple_merge[simple_merge$age > -1510000 & simple_merge$age <= -1492000, "MIS_stage"] <- "MIS_050"
    simple_merge[simple_merge$age > -1530000 & simple_merge$age <= -1510000, "MIS_stage"] <- "MIS_051"
    simple_merge[simple_merge$age > -1547500 & simple_merge$age <= -1530000, "MIS_stage"] <- "MIS_052"
    simple_merge[simple_merge$age > -1570000 & simple_merge$age <= -1547500, "MIS_stage"] <- "MIS_053"
    simple_merge[simple_merge$age > -1585000 & simple_merge$age <= -1570000, "MIS_stage"] <- "MIS_054"
    simple_merge[simple_merge$age > -1608000 & simple_merge$age <= -1585000, "MIS_stage"] <- "MIS_055"
    simple_merge[simple_merge$age > -1628500 & simple_merge$age <= -1608000, "MIS_stage"] <- "MIS_056"
    simple_merge[simple_merge$age > -1642500 & simple_merge$age <= -1628500, "MIS_stage"] <- "MIS_057"
    simple_merge[simple_merge$age > -1670000 & simple_merge$age <= -1642500, "MIS_stage"] <- "MIS_058"
    simple_merge[simple_merge$age > -1697500 & simple_merge$age <= -1670000, "MIS_stage"] <- "MIS_059"
    simple_merge[simple_merge$age > -1715000 & simple_merge$age <= -1697500, "MIS_stage"] <- "MIS_060"
    simple_merge[simple_merge$age > -1743000 & simple_merge$age <= -1715000, "MIS_stage"] <- "MIS_061"
    simple_merge[simple_merge$age > -1758000 & simple_merge$age <= -1743000, "MIS_stage"] <- "MIS_062"
    simple_merge[simple_merge$age > -1782000 & simple_merge$age <= -1758000, "MIS_stage"] <- "MIS_063"
    simple_merge[simple_merge$age > -1802500 & simple_merge$age <= -1782000, "MIS_stage"] <- "MIS_064"
    simple_merge[simple_merge$age > -1816000 & simple_merge$age <= -1802500, "MIS_stage"] <- "MIS_065"
    simple_merge[simple_merge$age > -1826000 & simple_merge$age <= -1816000, "MIS_stage"] <- "MIS_066"
    simple_merge[simple_merge$age > -1832500 & simple_merge$age <= -1826000, "MIS_stage"] <- "MIS_067"
    simple_merge[simple_merge$age > -1849000 & simple_merge$age <= -1832500, "MIS_stage"] <- "MIS_068"
    simple_merge[simple_merge$age > -1859500 & simple_merge$age <= -1849000, "MIS_stage"] <- "MIS_069"
    simple_merge[simple_merge$age > -1875000 & simple_merge$age <= -1859500, "MIS_stage"] <- "MIS_070"
    simple_merge[simple_merge$age > -1898000 & simple_merge$age <= -1875000, "MIS_stage"] <- "MIS_071"
    simple_merge[simple_merge$age > -1915000 & simple_merge$age <= -1898000, "MIS_stage"] <- "MIS_072"
    simple_merge[simple_merge$age > -1941000 & simple_merge$age <= -1915000, "MIS_stage"] <- "MIS_073"
    simple_merge[simple_merge$age > -1965000 & simple_merge$age <= -1941000, "MIS_stage"] <- "MIS_074"
    simple_merge[simple_merge$age > -1990000 & simple_merge$age <= -1965000, "MIS_stage"] <- "MIS_075"
    simple_merge[simple_merge$age > -2017000 & simple_merge$age <= -1990000, "MIS_stage"] <- "MIS_076"
    simple_merge[simple_merge$age > -2043000 & simple_merge$age <= -2017000, "MIS_stage"] <- "MIS_077"
    simple_merge[simple_merge$age > -2088000 & simple_merge$age <= -2043000, "MIS_stage"] <- "MIS_078"
    simple_merge[simple_merge$age > -2103000 & simple_merge$age <= -2088000, "MIS_stage"] <- "MIS_079"
    simple_merge[simple_merge$age > -2125000 & simple_merge$age <= -2103000, "MIS_stage"] <- "MIS_080"
    simple_merge[simple_merge$age > -2146000 & simple_merge$age <= -2125000, "MIS_stage"] <- "MIS_081"
    simple_merge[simple_merge$age > -2168000 & simple_merge$age <= -2146000, "MIS_stage"] <- "MIS_082"
    simple_merge[simple_merge$age > -2192000 & simple_merge$age <= -2168000, "MIS_stage"] <- "MIS_083"
    simple_merge[simple_merge$age > -2207500 & simple_merge$age <= -2192000, "MIS_stage"] <- "MIS_084"
    simple_merge[simple_merge$age > -2236000 & simple_merge$age <= -2207500, "MIS_stage"] <- "MIS_085"
    simple_merge[simple_merge$age > -2250000 & simple_merge$age <= -2236000, "MIS_stage"] <- "MIS_086"
    simple_merge[simple_merge$age > -2273000 & simple_merge$age <= -2250000, "MIS_stage"] <- "MIS_087"
    simple_merge[simple_merge$age > -2291000 & simple_merge$age <= -2273000, "MIS_stage"] <- "MIS_088"
    simple_merge[simple_merge$age > -2309000 & simple_merge$age <= -2291000, "MIS_stage"] <- "MIS_089"
    simple_merge[simple_merge$age > -2333000 & simple_merge$age <= -2309000, "MIS_stage"] <- "MIS_090"
    simple_merge[simple_merge$age > -2350000 & simple_merge$age <= -2333000, "MIS_stage"] <- "MIS_091"
    simple_merge[simple_merge$age > -2373000 & simple_merge$age <= -2350000, "MIS_stage"] <- "MIS_092"
    simple_merge[simple_merge$age > -2387000 & simple_merge$age <= -2373000, "MIS_stage"] <- "MIS_093"
    simple_merge[simple_merge$age > -2407000 & simple_merge$age <= -2387000, "MIS_stage"] <- "MIS_094"
    simple_merge[simple_merge$age > -2427000 & simple_merge$age <= -2407000, "MIS_stage"] <- "MIS_095"
    simple_merge[simple_merge$age > -2452000 & simple_merge$age <= -2427000, "MIS_stage"] <- "MIS_096"
    simple_merge[simple_merge$age > -2477000 & simple_merge$age <= -2452000, "MIS_stage"] <- "MIS_097"
    simple_merge[simple_merge$age > -2494000 & simple_merge$age <= -2477000, "MIS_stage"] <- "MIS_098"
    simple_merge[simple_merge$age > -2510000 & simple_merge$age <= -2494000, "MIS_stage"] <- "MIS_099"
    simple_merge[simple_merge$age > -2540000 & simple_merge$age <= -2510000, "MIS_stage"] <- "MIS_100"
    simple_merge[simple_merge$age > -2554000 & simple_merge$age <= -2540000, "MIS_stage"] <- "MIS_101"
    simple_merge[simple_merge$age > -2575000 & simple_merge$age <= -2554000, "MIS_stage"] <- "MIS_102"
    simple_merge[simple_merge$age > -2595000 & simple_merge$age <= -2575000, "MIS_stage"] <- "MIS_103"
    simple_merge$age <- simple_merge$age/-1000
    simple_merge <- simple_merge[complete.cases(simple_merge),]
    simple_merge$MIS_stage <- as.character(simple_merge$MIS_stage)
    return(simple_merge)
}
