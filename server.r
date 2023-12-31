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
      
      selectInput("proxy_type", "Proxy", choices=c("All", proxy_types), selected="isotope")
      
  })
  
  output$season_type_ui <- renderUI({
      
      selectInput("season_type", "Season", choices=season_types, selected="annual")
      
  })
  
  
  output$unit_type_ui <- renderUI({
      
      selectInput("unit_type", "Units", choices=unit_types, selected="permil")
      
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
      if(length(data_list) > 1){
          if(input$scale==TRUE){
              simple_merge <- data.frame(age=unlist(sapply(data_list, function(x) x$age)), value=unlist(sapply(data_list, function(x) scaleTransform(x$paleoData_values))), doi=as.character(unlist(sapply(data_list, function(x) rep(x$pub1_doi, length(x$paleoData_values))))))
          } else if(input$scale==FALSE){
              simple_merge <- data.frame(age=unlist(sapply(data_list, function(x) x$age)), value=unlist(sapply(data_list, function(x) x$paleoData_values)), doi=as.character(unlist(sapply(data_list, function(x) rep(x$pub1_doi, length(x$paleoData_values))))))
          }
      } else if(length(data_list)==1){
          simple_merge <- data.frame(age=data_list[["age"]], value=data_list[["paleoData_values"]], doi=rep(data_list[["pub1_doi"]], length(data_list[["paleoData_values"]])))
      } else if(length(data_list)==0){
          NULL
      }
      
      
      simple_merge <- simple_merge[simple_merge$doi %in% zipsInBounds()$pub1_doi,]
      simple_merge <- simple_merge[complete.cases(simple_merge),]
      
      simple_merge
      
  })

	dataFramePrep <- reactive({
		simple_merge <- dataFrameRaw()

		if(input$moving_average>0){
			simple_merge <- simple_merge[,c("age", "value")]
			#simple_merge$age <- simple_merge$age*100
          	simple_merge_years <- data.frame(age=seq(round(min(simple_merge$age), 0), round(max(simple_merge$age), 0), 1))
			simple_merge <- data.table::as.data.table(simple_merge)
			simple_merge_years <- data.table::as.data.table(simple_merge_years)
			data.table::setkey(simple_merge, age)
			data.table::setkey(simple_merge_years, age)
			simple_merge <- as.data.frame(merge(simple_merge, simple_merge_years, by = "age", all=T, fill=T))
			simple_merge <- simple_merge %>%
    			mutate(value = na.approx(value))

			#simple_merge$age <- simple_merge$age/100
      	}
        
        simple_merge
	
	})


	dataFrame <- reactive({
		
		if(input$moving_average>1){
		simple_merge <- dataFramePrep()
  			simple_merge$value <- TTR::DEMA(simple_merge$value, input$moving_average)
		} else if(input$moving_average==0){
			simple_merge <- dataFrameRaw()
		}

		#simple_merge <- mis_assign(simple_merge)
        simple_merge

	})

    climate_ranges <- reactiveValues(x = NULL, y = NULL)

            observeEvent(input$plot_climate_dblclick, {
                brush <- input$plot_climate_brush
                if (!is.null(brush)) {
                    climate_ranges$x <- c(brush$xmin, brush$xmax)
                    
                } else {
                    climate_ranges$x <- NULL
                }
            })
            


  
  climatePlot <- reactive({
      
      simple_merge <- dataFrame()
      
      y_lab <- if(input$scale==TRUE){
          paste0("Scaled ", input$unit_type)
      } else if(input$scale==FALSE){
          input$unit_type
      }
      
      ggplot(simple_merge, aes(age, value)) +
      geom_line() + 
      scale_y_continuous(y_lab) +
      scale_x_continuous("cal years BP", labels=scales::comma, limits=input$age_range) +
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
        stroke=FALSE, fillOpacity=0.4, fillColor=pal(colorData), group="markers") %>%
      addLegend("bottomleft", pal=pal, values=colorData, title=colorBy,
        layerId="colorLegend")
  })
  

  showZipcodePopup <- function(zipcode, lat, lng) {
    selectedZip <- metadata[metadata$pub1_doi == zipcode,]
    
    timeseries_list <- TS[which(doi==zipcode)]
    
    
    
	timeseries_data <- data.frame(age=as.numeric(timeseries_list[[1]]$age), value=as.numeric(timeseries_list[[1]]$paleoData_values))
    
    y_lab <- timeseries_list[[1]]$paleoData_proxy
    
	timeseries_plot <- ggplot(timeseries_data, aes(age, value)) +
		geom_line() +  
      	scale_y_continuous(y_lab) +
      	scale_x_continuous("cal years BP", labels=scales::comma) +
      	theme_light()
	
    fldr <- tempfile()
    dir.create(fldr)
    print(fldr)
    ggsave(filename = paste(fldr, "test.png", sep = "/"), timeseries_plot)
    img_path <- gsub("//", "/", paste(fldr, "test.png", sep = "/"))
    print(img_path)
    #tst <- paste(readLines(paste(fldr, "test.png", sep = "/")), collapse = "")
    #print(tst)

    selectedZip <- selectedZip[1,]
    content <- as.character(tagList(
      tags$h4(as.character(selectedZip$siteName)),
      paste0("Min: ", round(selectedZip$minYear, 0)), tags$br(),
      paste0("Max: ", round(selectedZip$maxYear, 0)), tags$br(),
      paste0("Resolution: ", round(selectedZip$agesPerKyr*1000), " years"),
      paste0("Archive Type: ", selectedZip$archiveType), tags$br(),
      paste0("Proxy Record: ", selectedZip$proxy), tags$br(),
      paste0("Seasonailty: ", selectedZip$seasonGeneral), tags$br()
      #paste0(includeHTML(paste0(fldr, "/test.svg")))
      ))
      
      next_content <- sprintf(
      paste0('<img src="', paste(fldr, 'test.png', sep = '/'), '" />')
      ) %>% lapply(htmltools::HTML)
      
      
    leafletProxy("map") %>% addPopups(lng, lat, content, layerId = zipcode)
    
    #leafpop::addPopupImages(paste0("<img src='", img_path, "' />"), group="markers")
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
  filename = function() { paste("northAmerica_subset", '.csv') },
  content = function(file
  ) {
      write.csv(reactiveZip(), file)
  }
  )
  
  
  
}
