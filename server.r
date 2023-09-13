library(leaflet)
library(RColorBrewer)
library(scales)
library(lattice)
library(dplyr)
library(ggplot2)
library(data.table)
library(DT)

function(input, output, session) {
    
    
    # Leaflet bindings are a bit slow; for now we'll just sample to compensate
    set.seed(100)
    
    
    zipdataInput <- reactive({
    metadata
    
    })
    
    
    # By ordering by centile, we ensure that the (comparatively rare) SuperZIPs
    # will be drawn last and thus be easier to see
    zipdataToo <- reactive(
    {zipdata <- zipdataInput()
    #zipdata[order(zipdata$centile),]
    })
    


  ## Interactive Map ###########################################

  # Create the map
  output$map <- renderLeaflet({
      leaflet() %>%
            addTiles() %>%
      setView(lng = -119.2, lat = 45.7, zoom = 3)
  })

  # A reactive expression that returns the set of zips that are
  # in bounds right now
  zipsInBounds <- reactive({
      zipdata <- zipdataToo()
    if (is.null(input$map_bounds))
      return(zipdata[FALSE,])
    bounds <- input$map_bounds
    latRng <- range(bounds$north, bounds$south)
    lngRng <- range(bounds$east, bounds$west)

    subset(zipdata,
      lat >= latRng[1] & lat <= latRng[2] &
        lon >= lngRng[1] & lon <= lngRng[2])
  })

  # Precalculate the breaks we'll need for the two histograms
  #centileBreaks <- hist(plot = FALSE, allzips$centile, breaks = 20)$breaks

  output$histCentile <- renderPlot({
    # If no zipcodes are in view, don't plot
    if (nrow(zipsInBounds()) == 0)
      return(NULL)

    hist(zipsInBounds()$centile,
      breaks = 20, #centileBreaks,
      main = "Chaco Isotope Data",
      xlab = "Percentile",
      #xlim = range(allzips$centile),
      col = '#00DD00',
      border = 'white')
  })
  
  output$archive_type_ui <- renderUI({
      
      selectInput("archive_type", "Archive", choices=c("All", archive_types), selected="All")
      
  })
  
  output$proxy_type_ui <- renderUI({
      
      selectInput("proxy_type", "Proxy", choices=c("All", proxy_types), selected="All")
      
  })
  
  output$season_type_ui <- renderUI({
      
      selectInput("season_type", "Season", choices=season_types, selected="annual")
      
  })
  
  
  output$unit_type_ui <- renderUI({
      
      selectInput("unit_type", "Units", choices=unit_types, selected="degC")
      
  })
  
  dataPull <- reactive({
      

      data_unit = which(ug == input$unit_type)
      data_season = which(sg == input$season_type)
      data_sub = data_unit[data_unit %in% data_season]

      if(input$proxy_type != "All"){
          data_proxy = which(pg == input$proxy_type)
          data_sub = data_unit[data_sub %in% data_proxy]
      }
      
      if(input$archive_type != "All"){
          data_archive = which(ag == input$archive_type)
          data_sub = data_unit[data_sub %in% data_archive]
      }
      
      
      TS[data_sub]

  })
  
  dataFrameRaw <- reactive({
      
      data_list <- dataPull()
      
      if(input$scale==TRUE){
          simple_merge <- data.frame(age=unlist(sapply(data_list, function(x) x$age)), value=unlist(sapply(data_list, function(x) scaleTransform(x$paleoData_values))), doi=unlist(sapply(data_list, function(x) rep(x$pub1_doi, length(x$paleoData_values)))))
      } else if(input$scale==FALSE){
          simple_merge <- data.frame(age=unlist(sapply(data_list, function(x) x$age)), value=unlist(sapply(data_list, function(x) x$paleoData_values)), doi=unlist(sapply(data_list, function(x) rep(x$pub1_doi, length(x$paleoData_values)))))
      }
      
      simple_merge <- simple_merge[simple_merge$doi %in% zipsInBounds()$pub1_doi,]
      simple_merge <- simple_merge[complete.cases(simple_merge),]
      
      simple_merge
      
  })

	dataFramePrep <- reactive({
		simple_merge <- dataFrameRaw()

		if(input$moving_average>0){
			simple_merge <- simple_merge[,c("age", "value")]
			simple_merge$age <- simple_merge$age*100
          	simple_merge_years <- data.frame(age=seq(round(min(simple_merge$age), 0), round(max(simple_merge$age), 0), 1))
			simple_merge <- data.table::as.data.table(simple_merge)
			simple_merge_years <- data.table::as.data.table(simple_merge_years)
			data.table::setkey(simple_merge, age)
			data.table::setkey(simple_merge_years, age)
			simple_merge <- as.data.frame(data.table::fset(simple_merge, simple_merge_years, by = "id"))
			simple_merge <- simple_merge %>%
    			mutate(value = na.approx(value))
      	}
	
	})


	dataFrame <- reactive({
		
		if(input$moving_average>1){
		simple_merge <- dataFramePrep()
  			simple_merge$value <- TTR::DEMA(simple_merge$value, input$moving_average)
		} else if(input$moving_average==0){
			simple_merge <- dataFrameRaw()
		}

		simple_merge <- mis_assign(simple_merge)

	})


  
  climatePlot <- reactive({
      
      simple_merge <- dataFrame()
      
      ggplot(simple_merge, aes(age, value)) +
      geom_line() + 
      scale_y_continuous(input$unit_type) +
      scale_x_continuous("kiloyears") +
      theme_light()
      
      
  })
  
  output$climate_plot <- renderPlot({
      climatePlot()
  })
  
  
  isoDensData <- reactive({
      
      data <- zipsInBounds()
      data_melt <- reshape2::melt(data[,c("siteName", "minYear", "maxYear")])
      
      dist.plot <- ggplot(data_melt, aes(value, colour=variable, fill=variable)) +
          geom_density(alpha=0.6) +
          theme_light() +
          facet_wrap(.~variable, scales="free", ncol=1) + 
          theme(legend.position="none")
      
      

      
      dist.plot
      
  })
  
  output$isoDens <- renderPlot({
      
      isoDensData()
      
  })

  output$scatterCollegeIncome <- renderPlot({
    # If no zipcodes are in view, don't plot
    if (nrow(zipsInBounds()) == 0)
      return(NULL)
      
      colorBy <- input$color
      sizeBy <- input$size
      
      cor.test <- lm(as.numeric(as.vector(zipsInBounds()[[colorBy]])) ~ as.numeric(as.vector(zipsInBounds()[[sizeBy]])))
      
      r2 <- summary(cor.test)$r.squared
      intercept <- cor.test$coef[1]
      slope <- cor.test$coef[2]
      

#print(xyplot(zipsInBounds()[[colorBy]] ~ zipsInBounds()[[sizeBy]], xlab=sizeBy, ylab=colorBy, type=c("p", "r"), main=expression(paste("y ="*paste(slope)*"x + "*paste(intercept)*", r"^2*paste(r2)))))

    temp.frame <- data.frame(as.numeric(as.vector(zipsInBounds()[[sizeBy]])), as.numeric(as.vector(zipsInBounds()[[colorBy]])))
    colnames(temp.frame) <- c("x", "y")
    temp.frame <- na.omit(temp.frame)
    
    scatter <- ggplot(aes(x, y), data=temp.frame) +
    geom_point(colour="blue") +
    stat_smooth(method="lm") +
    theme_light() +
    scale_x_continuous(sizeBy) +
    scale_y_continuous(colorBy) +
     annotate("text", label=lm_eqn(cor.test), x=0, y=Inf, hjust=0, vjust=1, parse=TRUE) 
    
    scatter
  })
  


  # This observer is responsible for maintaining the circles and legend,
  # according to the variables the user has chosen to map to color and size.
  observe({
     zipdata <- zipdataToo()
    colorBy <- input$color
    sizeBy <- input$size

    if (colorBy == "superzip") {
      # Color and palette are treated specially in the "superzip" case, because
      # the values are categorical instead of continuous.
      colorData <- ifelse(zipdata$centile >= (100 - input$threshold), "yes", "no")
      pal <- colorFactor("Spectral", colorData)
    } else {
      colorData <- zipdata[[colorBy]]
      pal <- colorFactor(
        palette = 'Dark2',
        domain = zipdata[,colorBy]
      )

    }

    #if (sizeBy == "superzip") {
    #  # Radius is treated specially in the "superzip" case.
    #  radius <- ifelse(zipdata$centile >= (100 - input$threshold), 30000, 3000)
    #} else {
    #  radius <- zipdata[[sizeBy]] / max(zipdata[[sizeBy]]) * 30000
    #}

    leafletProxy("map", data = zipdata) %>%
      clearShapes() %>%
      addCircleMarkers(~lon, ~lat, radius=20, layerId=~pub1_doi,
        stroke=FALSE, fillOpacity=0.4, fillColor=pal(colorData)) %>%
      addLegend("bottomleft", pal=pal, values=colorData, title=colorBy,
        layerId="colorLegend")
  })
  

  showZipcodePopup <- function(zipcode, lat, lng) {
    selectedZip <- metadata[metadata$pub1_doi == zipcode,]
    selectedZip <- selectedZip[1,]
    content <- as.character(tagList(
      tags$h4(as.character(selectedZip$siteName)),
      paste0("Min: ", round(selectedZip$minYear, 0)), tags$br(),
      paste0("Max: ", round(selectedZip$maxYear, 0)), tags$br(),
      paste0("Resolution: ", round(selectedZip$agesPerKyr*1000), " years"),
      tags$hr(),
      paste0("Archive Type: ", selectedZip$archiveType), tags$br(),
      paste0("Proxy Record: ", selectedZip$proxy)), tags$br(),
      paste0("Seasonailty: ", selectedZip$seasonGeneral)
      )
    leafletProxy("map") %>% addPopups(lng, lat, content, layerId = zipcode)
  }

  # When map is clicked, show a popup with city info
  observe({
    leafletProxy("map") %>% clearPopups()
    event <- input$map_marker_click
    if (is.null(event))
      return()

    isolate({
      showZipcodePopup(event$id, event$lat, event$lng)
    })
  })
  
  ####When zipcode is selected, show popup with city info
  
  #observe({
  #    leafletProxy("map") %>% clearPopups()
  #    event <- as.numeric(paste(input$yourzipcode))
  #    zipframe <- subset(zipcodes, zipcodes$zip_code==event)
      
      
      
  #    if (is.null(event))
  #    return()
      
  #    isolate({
  #        showZipcodePopup(event, zipframe$latitude, zipframe$longitude)
  #    })
  #})



  ## Data Explorer ###########################################
  
  reactiveZip <- reactive({
      
      smalls <- zipsInBounds()
      
      smalls
      
  })



  output$ziptable <- DT::renderDataTable({
      
      
    df <- reactiveZip()
    
    df <- df[,c("siteName", "archiveType", "proxy", "variableName", "units", "seasonGeneral", "minYear", "maxYear", "agesPerKyr", "pub1_doi")]

    DT::datatable(df)
  })
  
  
  output$downloaddata <- downloadHandler(
  filename = function() { paste("chacoData", '.csv', sep=',') },
  content = function(file
  ) {
      write.csv(reactiveZip(), file)
  }
  )
  
  
  
}
