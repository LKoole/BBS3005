# Install RCy3
if(!"RCy3" %in% installed.packages()){
  install.packages("BiocManager")
  BiocManager::install("RCy3")
}
library(RCy3)

#list of cytoscape apps to install
installation_responses <- c()

#list of app to install
cyto_app_toinstall <- c("clustermaker2", "enrichmentmap", "autoannotate", "wordcloud", "stringapp", "aMatReader")

cytoscape_version <- unlist(strsplit( cytoscapeVersionInfo()['cytoscapeVersion'],split = "\\."))
if(length(cytoscape_version) == 3 && as.numeric(cytoscape_version[1]>=3) 
   && as.numeric(cytoscape_version[2]>=7)){
  for(i in 1:length(cyto_app_toinstall)){
    #check to see if the app is installed.  Only install it if it hasn't been installed
    if(!grep(commandsGET(paste("apps status app=\"", cyto_app_toinstall[i],"\"", sep="")), 
             pattern = "status: Installed")){
      installation_response <-commandsGET(paste("apps install app=\"", 
                                                cyto_app_toinstall[i],"\"", sep=""))
      installation_responses <- c(installation_responses,installation_response)
    } else{
      installation_responses <- c(installation_responses,"already installed")
    }
  }
  installation_summary <- data.frame(name = cyto_app_toinstall, 
                                     status = installation_responses)
  
  knitr::kable(list(installation_summary),
               booktabs = TRUE, caption = 'A Summary of automated app installation'
  )
}

# Make sure that you are connected to Cytoscape 
cytoscapePing()

# Confirm that Cytoscape is installed and opened
cytoscapeVersionInfo ()

# To see all functions in RCy3
help(package=RCy3)

# CyREST API list all the functions available in a base distribution of cytoscape
cyrestAPI()  # CyREST API

# To see available commands
commandsAPI()  # Commands API

commandsHelp("help")

# Get image of network
setwd("/Users/lisakoole/Desktop/CytoscapeR")
filename <- file.path(getwd(),"230_Isserlin_RCy3_intro", "images","initial_example_network.png")
exportImage(filename, type = "pdf")

# query the String Database to get all interactions found for our set of top genes
commandsHelp("help string")
commandsHelp("help string protein query")

mergenetworks(c("Network 1", "Network 2"), "Merged Network",  
              nodeKeys=c("HGNC","query term"))
