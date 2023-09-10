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

load("data/wNAm1_0_0.RData")
metadata <- read.csv("data/wNAm_v1.csv")
metadata$ID <- formatC(make.names(metadata$siteName), width=5, format="d", flag="0")
archive_types <- unique(sapply(TS, function(x) x$archiveType))
proxy_types <- unique(unlist(sapply(TS, function(x) x$paleoData_proxyGeneral)))
site_names <- unique(lipdR::pullTsVariable(sTS,"geo_siteName"))
season_types <-  unique(lipdR::pullTsVariable(sTS,"interpretation1_seasonalityGeneral"))

varsColor <- c("archiveType", "proxy", "variableName", "seasonGeneral")




