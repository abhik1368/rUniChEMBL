library(RCurl)
library(jsonlite)

#' get.scid.sid
#' @title Get the source compound ids from another source compound id
#' @description a list of all src_compound_ids from all sources which are
#' CURRENTLY assigned to the same structure as a currently assigned 
#' query src_compound_id. 
#' The output will include query src_compound_id if it is a valid src_compound_id 
#' with a current assignment.
#' @name get.scid.sid 
#' @docType package
#' @param x : Input string Source compound id
#' @param y : Input integer Source id
#' @export
#' @examples
#' \donttest{
#' # Get source compound ids and source information
#' # Using ChEMBL ID and source
#' get.scid.sid("CHEMBL12",1)
#' # Using drugbank id and source
#' get.scid.sid("DB00789",2)
#' }
get.scid.sid <- function(x,y) {
    url <- sprintf("https://www.ebi.ac.uk/unichem/rest/src_compound_id/%s/%d",x,y)
    h <- getCurlHandle()
    d <- getURL(url, curl=h)
    status <- getCurlInfo(h)$response.code
    rm(h)
    if (status == 200) {
        data(SourceNames)
        sn<-SourceNames
        data<-fromJSON(d)
        df<-merge(data,sn,by.x="src_id",by.y="src_id")
        return(df)
    } else {
        return(NULL)
    }
}
#' get.sAll.sid
#' @title Get the all source compound ids from another source compound id
#' @description Obtain a list of all src_compound_ids from all sources 
#' (including BOTH current AND obsolete assignments) to the same structure 
#' as a currently assigned query 
#' src_compound_id.
#' @name get.sAll.sid 
#' @docType package
#' @param x : Input chemblid
#' @param y : Input integer source id
#' @export
#' @examples
#' \donttest{
#' # Get all source ids using ChEMBL id and source
#' get.sAll.sid("CHEMBL12",1)
#' # Using drugbank id and source
#' get.sAll.sid("DB00789",2)
#' }
get.sAll.sid <- function(x,y) {
    url <- sprintf("https://www.ebi.ac.uk/unichem/rest/src_compound_id_all/%s/%d",x,y)
    h <- getCurlHandle()
    d <- getURL(url, curl=h)
    status <- getCurlInfo(h)$response.code
    rm(h)
    if (status == 200) {
        data(SourceNames)
        sn<-SourceNames
        data<-fromJSON(d)
        df<-merge(data,sn,by.x="src_id",by.y="src_id")
        return(df)
    }else {
        return(NULL)
    }
}

#' get.mapping.full
#' @title Get full mapping between two sources.
#' @description Obtain a full mapping between two sources. Uses only currently 
#' assigned src_compound_ids from both sources.
#' @name get.mapping.full
#' @docType package
#' @param x : Input integer source id
#' @param y : Input integer source id
#' @export
#' @examples
#' \donttest{
#' # Get full mapping of PDBe and ChEMBL
#' get.mapping.full(3,1)
#' # Get full mapping of ZINC and ChEMBL
#' get.mapping.full(9,1) 
#' }

get.mapping.full <-function(x,y) {
    url <- sprintf("https://www.ebi.ac.uk/unichem/rest/mapping/%d/%d",x,y)
    h <- getCurlHandle()
    d <- getURL(url, curl=h)
    status <- getCurlInfo(h)$response.code
    rm(h)
    if (status == 200) {
        return(fromJSON(d))
    } else {
        return(NULL)
    }
}
#' get.src_id.InCHIKey
#' @title Get source compound ids
#' @description Obtain a list of src_compound_ids (from all sources) which 
#' are CURRENTLY  assigned to a query InChI Key. Returns a list of data from 
#' Unichem and ChEMBL databases.
#' @name get.src_id.InCHIKey 
#' @docType package
#' @param x : Input string InCHI Key
#' @export
#' @examples
#' \donttest{
#' # Get source compound ids from InCHIKey 
#' get.sid.InCHIKey("AAOVKJBEBIDNHE-UHFFFAOYSA-N")
#'
#' data<-get.sid.InCHIKey("BSYNRYMUTXBXSQ-UHFFFAOYSA-N")
#' # to get chembl data 
#' data$Chem
#' to get Unichem data
#' data$Uni
#' }
get.sid.InCHIKey<-function(x){
    url_uni <- sprintf("https://www.ebi.ac.uk/unichem/rest/inchikey/%s",x)
    url_chem<-sprintf("https://www.ebi.ac.uk/chemblws/compounds/stdinchikey/%s.json",x)
    h <- getCurlHandle()
    d <- getURL(url_uni, curl=h)
    c<-getURL(url_chem, curl=h)
    status <- getCurlInfo(h)$response.code
    rm(h)
    if (status == 200) {
        data(SourceNames)
        sn<-SourceNames
        data<-fromJSON(d)
        df<-merge(data,sn,by.x="src_id",by.y="src_id")
        cf<-do.call(rbind, lapply(fromJSON(c), data.frame))
        return(list(uni=df,Chem=cf))
    } else {
        return(NULL)
    }
}

#' get.sAll.InCHIKey
#' @title Get all src_compound_ids.
#' @description Get a list of all src_compound_ids (from all sources) which 
#' have current AND obsolete assignments to a query InChIKey
#' @name get.sAll.InCHIKey
#' @docType package
#' @param x : Input string InCHI Key
#' @export
#' @examples
#' \donttest{
#' # Get all the IDs using InCHIKey
#'  get.sAll.InCHIKey("AAOVKJBEBIDNHE-UHFFFAOYSA-N")
#' }

get.sAll.InCHIKey<-function(x){
    url <- sprintf("https://www.ebi.ac.uk/unichem/rest/inchikey_all/%s",x)
    h <- getCurlHandle()
    d <- getURL(url, curl=h)
    status <- getCurlInfo(h)$response.code
    rm(h)
    if (status == 200) {
        data(SourceNames)
        sn<-SourceNames
        data<-fromJSON(d)
        df<-merge(data,sn,by.x="src_id",by.y="src_id")
        return(df)
    }else {
        return(NULL)
    }
}

#' get.structure
#' @title  Get structure
#' @description Get structure(s) currently assigned to a query src_compound_id.
#' @name get.structure 
#' @docType package
#' @param x : Input string chemblid
#' @param s : Input integer source id (default is 1)
#' @export
#' @examples
#' \donttest{
#' # Get Standard inhci and InCHIKey from drugbank compound and source
#' get.structure("DB00321",2)
#' # Using ChEMBL compound and source id
#' get.structure("CHEMBL1231",1)
#' }

get.structure<-function(x,s=1){
    url <- sprintf("https://www.ebi.ac.uk/unichem/rest/structure/%s/%d",x,s)
    h <- getCurlHandle()
    d <- getURL(url, curl=h)
    status <- getCurlInfo(h)$response.code
    rm(h)
    if (status == 200) {
        return(fromJSON(d))
    } else {
        return(NULL)
    }
}

#' get.struc.all
#' @title Get Structures for source compound id
#' @description Get the standard InCHI and standard InCHI Key for the source 
#' compound id
#' @name get.struc.all
#' @docType package
#' @param x : Input string chemblid
#' @param s : Input integer source id (default is 1)
#' @export
#' @examples
#' \donttest{
#' # Get all the structure information using ChEMBL id and source.
#' get.struc.all("CHEMBL1231",1)
#' #using drugbank id and source
#' get.structure("DB00321",2)
#' }
get.struc.all<-function(x,s=1){
    url <- sprintf("https://www.ebi.ac.uk/unichem/rest/structure_all/%s/%d",x,s)
    h <- getCurlHandle()
    d <- getURL(url, curl=h)
    status <- getCurlInfo(h)$response.code
    rm(h)
    if (status == 200) {
        return(fromJSON(d))
    } else {
        return(NULL)
    }
}

#' get.url.sid
#' @title Get url for the query compound
#' @description Get a list of URLs for all src_compound_ids, from a specified 
#' source .
#' @name get.url.sid 
#' @docType package
#' @param x : Input string source compound id
#' @param y : Input integer source id
#' @param z : Input integer to source id
#' @export
#' @examples
#' \donttest{
#' # get urls of compounds using source compound id, source id
#' # get drugbank url from ChEMBL source id and ChEMBL source
#' get.url.sid("ChEMBL490",1,2)
#' 
#' # get chembl url from drugbank id and source 
#' get.url.sid("DB00715",2,1)
#' }
get.url.sid<-function(x,y,z){
    url <- sprintf("https://www.ebi.ac.uk/unichem/rest/src_compound_id_url/%s/%d/%d",x,y,z)
    h <- getCurlHandle()
    d <- getURL(url, curl=h)
    status <- getCurlInfo(h)$response.code
    rm(h)
    if (status == 200) {
        return(fromJSON(d))
    } else {
        return(NULL)
    }
}
#' get.sAll.obs
#' @title Get source compound id from obsolete source compound id
#' @description Get a list of all src_compound_ids from all sources with BOTH 
#' current AND obsolete to the same structure with an obsolete assignment to the #' query src_compound_id.
#' @name get.SrcAll.obs 
#' @docType package
#' @param x : Input string source compound id
#' @param y : Input integer to source id
#' @export
#' @examples
#' \donttest{
#' #get for drugbank compound and source 
#' get.sAll.obs("DB07699",2)
#' #get for chembl compound and source
#' get.sAll.obs("CHEMBL12",1)
#' }
get.sAll.obs<-function(x,y){
    url <- sprintf("https://www.ebi.ac.uk/unichem/rest/src_compound_id_all_obsolete/%s/%d",x,y)
    h <- getCurlHandle()
    d <- getURL(url, curl=h)
    status <- getCurlInfo(h)$response.code
    rm(h)
    if (status == 200) {
        data(SourceNames)
        sn<-SourceNames
        data<-fromJSON(d)
        df<-merge(data,sn,by.x="src_id",by.y="src_id")
        return(df)
    } else {
        return(NULL)
    }
}
#' get.verbose.InCHIkey
#' @title Get all src_compound_ids to a query InChIKey
#' @description Returns a dataframe containing src_id (the src_id for this source), 
#' src_url (the main home page of the source), name (the unique name for the source in 
#' UniChem, always lower case), name_long (the full name of the source, as defined by the 
#' source),name_label (A name for the source 
#' suitable for use as a 'label' for the source within a web-page. Correct case setting 
#' for source, and always less than 30 characters), description (a description of the 
#' content of the source), base_id_url_available (an flag indicating whether this source 
#' provides a valid base_id_url for creating cpd-specific links [1=yes, 0=no]),base_id_url 
#' (the base url for constructing hyperlinks to this source [append an identifier from 
#' this source to the end of this url to create a valid url to a specific page for this 
#' cpd], unless aux_for_url=1), aux_for_url (A flag to indicate whether the aux_src field 
#' should be used to create hyperlinks instead of the src_compound_id [1=yes, 0=no] , 
#' src_compound_id (a list of src_compound_ids from this source which are currently 
#' assigned to the query InChIKey, aux_src (a list of src-compound_id keys mapping to 
#' corresponding auxiliary data (url_id:value), for creating links if aux_for_url=1. Only 
#' shown if aux_for_url=1).
#' @name get.verbose.InCHIkey
#' @docType package
#' @param x : Input string InCHI Key
#' @export
#' @examples
#' \donttest{
#' # get for InCHIkey 
#' get.verbose.InCHIkey("GZUITABIAKMVPG-UHFFFAOYSA-N") 
#' 
#' get.verbose.InCHIkey("AAOVKJBEBIDNHE-UHFFFAOYSA-N") 
#' }

get.verbose.InCHIkey<-function(x){
    url <- sprintf("https://www.ebi.ac.uk/unichem/rest/verbose_inchikey/%s",x)
    h <- getCurlHandle()
    d <- getURL(url, curl=h)
    status <- getCurlInfo(h)$response.code
    rm(h)
    if (status == 200) {
        return(fromJSON(d))
    } else {
        return(NULL)
    }
}
#' get.cmp.inf
#' @title Get compound information from ChEMBL
#' @description These functions allow one to retrieve compounds information from ChEMBL
#' compounds are identified either by a ChEMBL ID or by a standard InChI key.
#' @name get.cmp.inf
#' @docType package
#' @param x : String representing chemblid or standard InCHI key for the molecule.
#' @param type : For \code{get.compound}, one of \code{chemblid} or 
#' \code{stdinchi} to indicate the nature of the molecule id. 
#' @export
#' @examples
#' \donttest{
#' #get information for chembl compound id
#' get.compound("CHEMBL12")
#' 
#' #get information for standard inchi
#' get.compound("QFFGVLORLPOAEC-SNVBAGLBSA-N",type='stdinchi')
#' }
get.cmp.inf <- function(x, type='chemblid') {
    types <- c('chemblid', 'stdinchi')
    type <- pmatch(type, types)
    if (is.na(type)) stop("Invalid type given")
    url <- switch(type,
                  url = 'https://www.ebi.ac.uk/chemblws/compounds/',
                  url = 'https://www.ebi.ac.uk/chemblws/compounds/stdinchikey/')
    url <- sprintf('%s%s.json', url, id)
    h <- getCurlHandle()
    d <- getURL(url, curl=h)
    status <- getCurlInfo(h)$response.code
    rm(h)
    if (status == 200) {
        return(fromJSON(d))
    } else {
        return(NULL)
    }
}

#' get.cmp.sim
#' @title Retrive similar compounds from ChEMBL database.
#' @description This function retrieves a dataframe of similar compounds 
#' from ChEMBL database given a smiles string as query and also given a similarity score above 70.
#' @name get.cmp.sim 
#' @docType package
#' @param mol : String representing smiles of the moelcule
#' @param sim : Integer representing for percentage of similarity 
#' for the query compound and the database molecules. Values ranges 
#' from 70 to 100.
#' @export
get.cmp.sim <- function(mol,sim=70) {
    url <- 'https://www.ebi.ac.uk/chemblws/compounds/similarity/'
    url <- sprintf('%s%s/%d.json', url,mol,sim)
    h <- getCurlHandle()
    d <- getURL(url, curl=h)
    status <- getCurlInfo(h)$response.code
    rm(h)
    if (status == 200) {
        return(fromJSON(d)[[1]])
    } else {
        return(NULL)
    }  
}
#' get.cmp.substruct
#' @title Get compound information from substructure query smiles.
#' @description This function retrieves a dataframe of all compounds from ChEMBL database 
#' containing the substructure represented by the given Canonical SMILES and their  
#' chemical properties.
#' @name get.compound.substruct 
#' @docType package
#' @param mol : String representing smiles of the moelcule
#' @export
#' @examples
#' \donttest{
#' #get compounds by substructure
#' get.cmp.subsruct("CN(CCCN)c1cccc2ccccc12")
#' }

get.cmp.substruct<-function(mol){
    url <- sprintf('https://www.ebi.ac.uk/chemblws/compounds/substructure/%s.json',mol)
    h <- getCurlHandle()
    d <- getURL(url, curl=h)
    status <- getCurlInfo(h)$response.code
    rm(h)
    if (status == 200) {
        data<-fromJSON(d)
        return(data$compounds)
    } else {
        return(NULL)
    }  
}
#' get.appDrugs
#' @title Get approved drugs for target.
#' @description  This function retrieves a dataframe of all approved drug compounds from
#' ChEMBL database given a string of ChEMBL target ID.
#' @name get.appDrugs
#' @docType package
#' @param x : string  ChEMBL target ID.
#' @export
#' @examples
#' \donttest{
#' #get chembl ids of approved drugs
#' get.appDrugs("CHEMBL1824")
#' }
get.appDrugs<-function(x){
    url<-sprintf('https://www.ebi.ac.uk/chemblws/targets/%s/approvedDrug.json',x)
    h <- getCurlHandle()
    d <- getURL(url, curl=h)
    status <- getCurlInfo(h)$response.code
    rm(h)
    if (status == 200) {
        return(fromJSON(d))[[1]]
    } else {
        return(NULL)
    }  
}
#' get.bioactivity.
#' @title Get Bioactivity Information for Compounds, Targets or Assays.
#' @description This method retrieves bioactivity information for a compound 
#' across multiple targets/assays or else for a target across multiple compounds. 
#' The function can also be used to retrieve all activities within a given assay. 
#' In all cases, ChEMBL identifiers must be used.
#' @name get.bioactivity
#' @docType package
#' @param x : Input string chemblid
#' @param type : Input string \code{'compound'},\code{'target'},\code{'assay'}. Default is 
#' \code{'compound'}.
#' @export
#' @examples
#' \donttest{
#' # get bioactivities of compounds
#' get.bioactivity("CHEMBL12",type='compound')
#' 
#' # get compound bioactivities for targets
#' get.bioactivity("CHEMBL240",type="target")
#' 
#' # get bioactivities by assay
#' get.bioactivity("CHEMBL1217643",type='assay')
#' }
get.bioactivity <- function(x, type='compound') {
    types <- c('compound', 'target', 'assay')
    type <- pmatch(type, types)
    if (is.na(type)) stop("Invalid type given")
    url <- switch(type,
                  url = 'https://www.ebi.ac.uk/chemblws/compounds/%s/bioactivities.json',
                  url = 'https://www.ebi.ac.uk/chemblws/targets/%s/bioactivities.json',
                  url = 'https://www.ebi.ac.uk/chemblws/assays/%s/bioactivities.json')
    url <- sprintf(url,x)
    h <- getCurlHandle()
    d <- getURL(url, curl=h)
    status <- getCurlInfo(h)$response.code
    rm(h)
    if (status == 200) {
        data<-do.call(rbind, lapply(fromJSON(d), data.frame))
        return(data)
    } else {
        return(NULL)
    }
}
#' get.moa
#' @title Get mechanism of action 
#' @description This function retrieves a data frame of compounds and its 
#' mode of action for a compound (where compound is a drug) and drug targets.
#' @name get.moa 
#' @docType package
#' @param x : Input string chemblid
#' @export
#' @examples
#' \donttest{
#' # get moa of drug
#' get.moa("CHEMBL1642")
#' } 
get.moa<-function(x){
    url <- 'https://www.ebi.ac.uk/chemblws/compounds/'
    url<-sprintf('%s%s/drugMechanism.json',url,x)
    h <- getCurlHandle()
    d <- getURL(url, curl=h)
    status <- getCurlInfo(h)$response.code
    rm(h)
    if (status == 200) {
        data<- do.call(rbind, lapply(fromJSON(d), data.frame))
        return(data)
    } else {
        return(NULL)
    }  
}
#' get.targets
#' @title Get target information.
#' @description This function retrieves the target information by chembl id and 
#' uniprot id and also retrieves all the names of targets by organisms. When 
#' \code{org="Homo sapiens"} subsets the data frame by organism homo sapiens and retrieves 
#' all the Homo sapiens taregts
#' @name get.targets  
#' @docType package
#' @param x : Input string chemblid
#' @param type : Input string 'chemblid' or 'uniprot'
#' @param org : Input string species name like "Homo sapiens","Plasmodium falciparum" and etc.
#' @export
#' @examples
#' \donttest{
#' #get target information by chembl ids
#' get.targets("CHEMBL1862",type='chemblid')
#' 
#' #get target information by uniprot ids
#' get.targets("Q13936",type='uniprot')
#' 
#' #get all the target information using organism name
#' get.targets(org="Homo Sapiens")
#' }


get.targets <- function(x,type='chemblid',org=NULL){
    types <- c('chemblid', 'uniprot')
    type <- pmatch(type, types)
    if (is.na(type)) stop("Invalid type given")
    if(is.null(org)){
        url<-switch(type,'https://www.ebi.ac.uk/chemblws/targets/%s.json',
                    'https://www.ebi.ac.uk/chemblws/targets/uniprot/%s.json')
        url <- sprintf(url,x)
        h <- getCurlHandle()
        d <- getURL(url, curl=h)
        status <- getCurlInfo(h)$response.code
        rm(h)
        if (status == 200) {
            data<- do.call(rbind, lapply(fromJSON(d), data.frame))
            return(data)
        } else {
            return(NULL)
        }
    }else{
        url<-'https://www.ebi.ac.uk/chemblws/targets.json'
        h <- getCurlHandle()
        d <- getURL(url, curl=h)
        status <- getCurlInfo(h)$response.code
        rm(h)
        if (status == 200) {
            df <- do.call(rbind, lapply(fromJSON(d), data.frame))
            data<-df[ which(df$organism==org)]
            return (data)
        } else {
            return(NULL)
        }
    } 
    
}

