

library(RMySQL)



aspDbFetch <- function(query = NULL, user=user, password=pw, dbname=db, host=host_ip, port = con_port){
    ################################################################################
    # Provide query and get your result as dataframe, closes connection afterwards #
    ################################################################################

    if(is.null(query)==FALSE){
        con <- dbConnect(MySQL(), user=user, password=password, dbname=dbname, host=host, port = port)
        #print(paste0("Connection established to ", dbname))
        aspData <- dbGetQuery(con,query)
        #print(sprintf("Fetched result with %s rows and %s columns", dim(aspData)[1],dim(aspData)[2]))
        dbDisconnect(con)
        #print("Closing connection")
        return(aspData)

    }else{
        print("Please check your query")
    }
}
