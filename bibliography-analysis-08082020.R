###########################################################
####
####  STEP 0: Prepare the R-Studio Environment
####
###########################################################

# SET THE WORKING DIRECTORY:
  # These instructions assume that all R scripts, data files, etc. are located within this project directory.
  # Set the working directory to your project directory.  
  # To do this: Session-> Set Working Directory -> To Project Directory

# INSTALL REQUIRED R PACKAGES
  # Be sure you have the following R packages installed.  If not install those you do not have.
    # install.packages("bibliometrix")
    # install.packages("dplyr")

# LOAD REQUIRED 
  #Load the following packages - need to do this every session
    library(bibliometrix)
    library(dplyr)


  
  ###########################################################
  ####
  ####  STEP 1: Create a List of Authors
  ####
  ###########################################################

  #*****
  # Work with your dean of faculty's office to obtain a list of faculty for the time period(s)
  # you want to analyze.
  # This will probably be a list of active faculty members.  It could be limited to specific 
  # Department and Program
  # groups. Starting with a complete list in one place is very helpful.
  #*****
  
  ### Example: Short list of Bates Faculty Members
        #- Ryan Bavis
        #- Geneva Laurita
        #- Holly Ewing


###########################################################
####
####  STEP 2: CREATE ".bib" FILES CONTAINING BIBLIOGRAPHIC INFORMATION OF INDIVIDUAL AUTHORS
####
###########################################################

#Note: Bib info can come from Scopus, PubMed, WOS, or ISI WoK databases.
#See: https://cran.r-project.org/web/packages/bibliometrix/vignettes/bibliometrix-vignette.html

#USING SCOPUS AS THE SEARCH TOOL:
  #Perform "Author" search for each faculty member of a given department.
  #Export citation information for all publications from Scopus in the BibTex format (i.e. .bib files).
  #Note: Scopus will only export up to 2000 files at a single time.
  #Copy all Scopus.bib export files (unique names for each) into the working directory.


  ###-Suggested .bib File Naming Convention:
    #lastname-firstinitial-department-DDMMYY.bib
    #e.g. laurita-g-bateschem-022619.bib
  
  ### Example: Short list of Bates Faculty Members
    #- Ryan Bavis: "bavis-r-030419.bib"
    #- Geneva Laurita: "laurita-g-022619.bib"
    #- Holly Ewing: "ewing-h-030419.bib"
  
  ### Transfer all of these .bib files to the working directory:
    # e.g. "G:/My Drive/R Project Files/Bibliometrix Journal Analysis"
  
  ### Combine the .bib file names into a single variable called "files" to 
    # be used below for reading in the .bib file data.
    # An example for the format of this file name string is given below
    # Note that each file name should be surrounded by quotation marks
    # the quoted file names are separated by commas,
    # and the entire series is surrounded by open and closed parentheses. 
  
  ### Example of Bates Authors .bib file name string:
   
   files <-c("bavis-r-030419.bib","laurita-g-022619.bib","ewing-h-030419.bib")
   
###########################################################
####
####  STEP 3: READ INDIVIDUAL OR GROUPS OF ".bib" FILES INTO A DATAFRAME "M"
####
###########################################################

   ###########################################################
   ####
   ####  STEP 3A: Read the .bib file(s) exported from Scopus author publication results into a dataframe.
   ####
   ###########################################################
   
   # Syntax for reading/combining multiple .bib files is:
     # M <- convert2df("file1.bib","file2.bib", ...)
     # So, you can copy and paste the Authors .bib "files" name string from Step 2 directly into
     # this convert2df function.
   
   # Note the the arguments for this function are dependent upon the original bibliographic data source
   # (i.e. Scopus, Web of Science, PubMed)
     # From the bibliometrix documentation (https://cran.r-project.org/web/packages/bibliometrix/bibliometrix.pdf):
       # convert2df creates a bibliographic data frame with cases corresponding
       # to manuscripts and variables to Field Tag in the original export file.
       # convert2df accepts two additional arguments: dbsource and format.
       # The argument dbsource indicates from which database the collection has been downloaded.
       # The argument format indicates the original exported file type (i.e. bibtex).
     # Also, see: https://cran.r-project.org/web/packages/bibliometrix/vignettes/Data-Importing-and-Converting.html
   
  
   #####
   ## MATHEMATICS AND NATURAL SCIENCES [Small Example Set]
   #####

  M <- convert2df(file = files, dbsource = "scopus", format = "bibtex")

  
  ## Use the following line for ISI Web of Science data:
  # M <-convert2df(D, dbsource = "isi",format = "bibtex") #Remove comment mark to run
  
  
  
  ###########################################################
  ####
  ####  STEP 3B: [Optional and Not Really Necessary] Export Dataframe "M" to a .csv File
  ####
  ###########################################################
  
  
  # Specific Example of function syntax:
    # Create an .CSV file named "Dataframe-M-08062020.csv" with the data from 
    # dataframe "M".  The names of each column will be included with these data.
  
  write.csv(M, "Dataframe-M-08062020.csv", row.names = FALSE)

 
  
  ###########################################################
  ####
  ####  STEP 3C: [Optional and Not Really Necessary] Create Archived "M Dataframe" For Each Analysis
  ####
  ###########################################################
  
  # If you intend to analyze several collections of .bib files, it can be useful to save a copy of the
    # M dataframe as a different dataframe at this point.
    # This allows you to pick up the analysis at this point for each of these dataframes and saves 
    # the work of having to generate them again.

   M.Example.080620.df <- M


  ###########################################################
  ####
  ####  STEP 3D: [OPTIONAL] CREATE A SUBSET OF DATAFRAME "M" SCOPED BY AFFILIATE PUBLICATION YEAR RANGE
  ####
  ###########################################################

  # DEFINE FUNCTION "LIMIT_PY" FOR SUBSETTING DATA BY PUBLICATION YEAR

    # LIMIT_PY is a function that scopes publication data from
    # dataframe M between two specified years.
    # By default publication date range is set to 1800 and 2100.
    # This should include everything published by the selected authors
    # in the parent database.
  
    # To select a different publication date range, specify the
    # start and end year as arguments for the function.
    # For example: >LIMIT_PY(2009,2019) will subset the data in M
    # to include only affiliate publications with publication dates
    # from 2009-2019, inclusive.
  
  LIMIT_PY <- function(minyear=1800, maxyear=2100) {
    M_SCOPE_PY <<- (M[M$PY>=minyear & M$PY<=maxyear,]) # "<<-" means assign("M_SCOPE_PY", M_SCOPE_PY, envir = .GlobalEnv)
    # str(M_SCOPE_PY) # Remove the comment symbol at the start of this line to report the M_SCOPE_PY structure
  }

  ### Call the "LIMIT_PY" function to create a dataframe called "M_SCOPE_PY"
  
    # "M_SCOPE_PY will contain a subset of "M" dataframe scoped by Author publication date
    # (i.e. only articles published by your faculty within the specified year range will be included)
    # Syntax for the function is: LIMIT_PY(FirstYear,LastYear)
    # where "FirstYear" and "LastYear" are the first and last years, respectively, of the desired date range.

  LIMIT_PY(2010,2019) 
 

  ### Create Archived "M_SCOPE_PY Dataframe" scoped by date range for each analysis
  
    # Because M_SCOPE_PY is a generic dataframe that will be overwritten every time you run the LIMIT_PY
    # function (e.g. if you are interested in analyzing more than one year range), you can "save" this dataframe
    # by creating another with a differnt name.
  
    # Suggested naming convention:
    # M.[descriptor].[MMDDYY].[FirstYear]_[LastYear].df
   
  M.Example.102519.2010_2019.df <- M_SCOPE_PY

  ### Redefine M as the new subsetted M_SCOPE_PY to avoid renameing M throughout below.
  
    # Subsequent analysis of the M-type dataframe data assumes this dataframe is called "M".
    # The simplest way to proceed to "STEP 4: BIBLIOMETRIC ANALYSIS OF BIBLIOGRAPHIC DATAFRAME "M"" 
    # from this point is to name any dataframe you will carry forward to"M".

  M <- M_SCOPE_PY 

###########################################################
####
####  STEP 4: CREATE A DATAFRAME CONTAINING THE NUMBER OF TIMES AFFILIATE AUTHORS HAVE "PUBLISHED IN" A JOURNAL
####
###########################################################

  ###########################################################
  ####
  ####  STEP 4A: CREATE "pub_in.df" DATAFRAME
  ####
  ###########################################################
  
  # CREATE A DATAFRAME WHERE ROWS ARE SUBSETTED BY "PUBLISHED IN" JOURNAL TITLE
    # NOTE: You can generate results using either the "SO" or "JI" columns of dataframe "M"
    # which contain either the full or abbreviated journal titles, respectively.

#library(dplyr) % This package was probably already loaded above. But if not, mount it now.
  
pub_in.df <- M %>%
  #group_by(JI) %>% # Gives ISO Source Abbreviation
  #group_by(SO) %>% #Gives full length title
  group_by(SO,SN) %>% #Gives full length title, ISSN
  #group_by(SO,JI,SN) %>% #Gives full length title, abbreviated title, ISSN
  tally() # Tallies all instances of a given title (i.e. "JI" or "SO" in column "JI" or "SO" of dataframe "M", respectively.)

  # ADD COLUMN NAMES TO THE "pub_in.df" DATAFRAME

colnames(pub_in.df) <- c("Published_In", "ISSN", "Frequency") 

  ###########################################################
  ####
  ####  STEP 4B: [OPTIONAL AND NOT NECESSARY] CREATE "pub_in.df" DATAFRAME
  ####
  ###########################################################

  # ORDER DATA BY DECREASING FREQUENCY COUNT
    # pub_in.df has all the data, but it is not ordered.
    # This orders pub_in.df by citation frequency in descending order and creates pub_in2.df

pub_in2.df <-pub_in.df[order(-pub_in.df$Frequency),]

  ###########################################################
  ####
  ####  STEP 4C: [Optional and Not Really Necessary] WRITE THE "pub_in2.df" DATA TO .CSV FILE
  ####
  ###########################################################

write.csv(pub_in2.df, file = "pubIn_xxfac_allTime_AllPublishers_111219.csv", row.names = FALSE)

  # [OPTIONAL] MODIFY ISSN VALUES IN DATAFRAME TO FORCE MS EXCEL TO READ ISSNs IN .CSV FILE AS CHARACTERS 
    # NOTE: If you will use Excel to work with the CSV file, you must force the ISSN
    # values to be read as characters (i.e. not numbers).
    # Make sure Excel does not treat the ISSN column values as numbers.  When this happens, 
    # if the ISSN has any 0's in the first set of the 8 postitions (e.g. 00XXXXXX) 
    # the 0's get removed (i.e. XXXXXX).
  
    # This surrounds the ISSN values with quotes and then 
    # prepends an equal sign on to that (e.g. 00XXXXXX becomes ="00XXXXXX".)  
    # Excel reads this in as a character.

pub_in2.df$ISSN = paste0('="', pub_in2.df$ISSN, '"')
write.csv(pub_in2.df, file = "pubIn_xxfac_allTime_AllPublishers_111219.csv", row.names = FALSE)



###########################################################
####
####  STEP 5: CREATE A DATAFRAME CONTAINING THE NUMBER OF TIMES AFFILIATE AUTHORS
####          HAVE CITED A JOURNAL (I.E. "CITED IN" BIBLIOGRAPHY)
####
###########################################################

  ### STEP 5A.1: CREATE NEW DATAFRAME "RefList" BY EXTRACTING SOURCE TITLES FROM CR (I.E. CITED REFERENCES) DATA IN "M" DATAFRAME
    # Extract a list of source titles from the cited reference, "CR_SO", variable in the "M" 
    # dataframe using the metaTagExtraction function (i.e. a list of the sources cited by 
    # articles in this local library)
  
    # The "metaTagExtraction" function is described on pages 39-40 of the 'bibliometix' users manual 
    # (Version 3.0.2; https://cran.r-project.org/web/packages/bibliometrix/bibliometrix.pdf).
    # The following function description comes from that document.
    
    #######################
    ## metaTagExtraction Description from the bibliometrix.pdf
    ## This function extracts other field tags, different from the standard ISI/SCOPUS codify.
    ##
    ## metaTagExtraction Usage
    ## metaTagExtraction(M, Field = "CR_AU", sep = ";", aff.disamb = TRUE)
    ##
    ## Arguments
    ### M is a data frame obtained by the converting function convert2df. It is a data matrix 
    ### with cases corresponding to articles and variables to Field Tag in the original ISI or SCOPUS file.
    ##
    ### Field is a character object. New tag extracted from aggregated data is specified 
    ### by this string. 
    ### Field can be equal to one of this tags:
    ##### "CR_AU" First Author of each cited reference
    ##### "CR_SO" Source of each cited reference
    ##### "AU_CO" Country of affiliation for each co-author
    ##### "AU1_CO" Country of affiliation for the first author
    ##### "AU_UN" University of affiliation for each co-author and the corresponding author (AU1_UN)
    ##### "SR" Short tag of the document (as used in reference lists) sep is the field separator character. 
    ####### This character separates strings in each column of the data frame. The default is sep = ";".
    ##
    ##### aff.disamb is a logical. If TRUE and Field="AU_UN", then a disambiguation algorithm is used 
    ####### to identify and match scientific affiliations (univ, research centers, etc.).
    ####### The default is aff.disamb=TRUE.
    #######################  
  
  
    ## Extract the source (e.g. journal title) "CR_SO" of every reference cited
    ## in the affiliate publications (i.e. "CR" variable defined in dataframe "M").
    ## These extracted cited reference source data are appended to dataframe "M" as variable CR_SO to create
    ## a new dataframe named "RefList".
  
    # As configured below, the metaTagExtraction function will pull a complete list of 
    # source titles of all the references cited in the scholarly output from affiliate author searches
    # used to create the initial .bib files (and, therefore, the "M" dataframe).

RefList <- metaTagExtraction(M, Field = "CR_SO", sep = ";")

 
  ### [OPTIONAL] Create Archived "RefList Dataframe" for each analysis
    # RefList is a generic dataframe that will be overwritten every time you run the metaTagExtraction
    # as describe above.  It can be useful to "save" the data in this specific dataframe
    # by saving it to another dataframe with another name.
  
    # Suggested naming convention:
    # RefList.[descriptor].[MMDDYY].df
   
    RefList.Example.102519.df <- RefList
    
    View(RefList.Example.102519.df)


  ### STEP 5A.2: CREATE SOURCE TITLE FREQUENCY TABLE "so" AND READ THAT TABLE INTO A NEW DATAFRAME "so.df"
    # Extract a list of source titles from the cited reference, "CR_SO", variable in the "M" 
    # dataframe using the metaTagExtraction function (i.e. a list of the sources cited by 
    # articles in this local library)


  # CREATE A FREQUENCY TABLE FROM THE CITED REFERENCE SOURCE VARIABLE IN THE "RefList" DATAFRAME.
    # The tableTag function allows you to create a 2 column (i.e. Tab, Freq) frequency table
    # containing the cited source titles and the number of times they are cited.
  
    # Note that it looks at the RefList CR_SO column/variable generated by the
    # megaextraction function.  During this extraction,
    # all of the source titles are separated by a semicolon.
    # The tableTag function recognizes that the titles are separated
    # by a semicolon with the "sep" argument.

so <- tableTag(RefList, "CR_SO", sep = ";") #Generates "so" freq table.


  # READ "so" SOURCE FREQUENCY DATA TABLE AS A DATAFRAME
    # Create a dataframe from the "so" table data. The "so" table is read 
    # in as the dataframe "so.df". 
  
so.df <- as.data.frame(so)
  
  # Redefine so.df dataframe column (i.e. variable) names "Tab" and "Freq" to be "RefSource"
  # and "Frequency", respectively.
  
colnames(so.df) <- c("RefSource", "Frequency")

  # VIEW THE so.df DATAFRAME. 
  # This should contain two columns of paired data: a journal/reference 
  # source title and the frequency that title shows up in the dataframe.  
    # Note: because there are many variations on journal title abbreviations, some journals will
    # be displayed several times.  For example, in this sample dataset Proceedings of the National Academy
    # of Sciences of the United States of America is listed several different ways (e.g. 
    # PROC NATL ACAD SCI U S A (20 times), PROC NATL ACAD SCI USA (14 times), PROC NATL ACAD SCI (6 times),
    # PROCEEDINGS OF THE NATIONAL ACADEMY OF SCIENCES USA (2 times)).  This is a problem.

View(so.df)

  ###########################################################
  ####
  ####  STEP 5B: ADD ISSN VALUES TO SOURCE FREQUENCY DATAFRAME "so.df" TO CREATE "title_issn_pairing.df" DATAFRAME
  ####
  ###########################################################

      ########    
      ## DESCRIPTION ON INTENT: COMBINE ALTERNATE VERSIONS OF THE SAME SOURCE TITLE
      ##   The Dream:
      ##   Use a "lookup table" containing all permutations of each journal title
      ##   and their corresponding ISSN to assign a "unique identifier" to
      ##   each source in the so.df frequency table.  This will allow us to assign
      ##   a single, unique identifier to all of the variations of a single journal
      ##   (e.g. Proceedings of the National Academy of Sciences and PNAS would both
      ##   be assigned the same ISSN, because they are the same journal.)
      ##   Of course, the real trick is building a lookup table with all (or enough) 
      ##   title-variant/ISSN pairs to be useful.
      ########

  # USE "addNewData" (and "readNewData") FUNCTION TO ADD ISSN VALUES FROM "title_issn_master.csv"
  # LOOKUP TABLE "so.df" TO GENERATE "title_issn_pairing.df"
      # NOTE: Unfortunately, because there are so many title/title-variants, it is very possible 
      #       that not all journal titles will be associated with an ISSN in the "title_issn_master.csv"
      #       lookup table. After running the addNewData/readNewData functions, identify 
      #       significant journal titles that are missing ISSN's and
      #       add the journal title/ISSN pair to the "title_issn_master.csv" lookup table.
      
      # NOTE: I created two additional R scripts to aid in the creation and
      #       maintenance of the lookup table.
      # NOTE: "scopusSourcesCsvToLookup_052319.R" converts the Scopus title list 
      #       file to a lookup table 
      # NOTE: "combine_csv_files_071019.R" to merge separate, lookup table 
      #       formated .csv files located in a subdirectory.
      #       This will allow you to combine lookup tables from other sources (friends!).
      
      # Note: The arguments for the addNewData function are ([name of .csv lookuptable],
      #       [title frequency dataframe],[the "allowedVars" variable]).
  
    ## DEFINE FUNCTIONS FOR LOOKUP TABLE ANALYSIS
      ### I got the original code (i.e. the functions "addNewData" and "readNewData" for 
      ### using a lookup table this way from:
      ### https://nicercode.github.io/blog/2013-07-09-modifying-data-with-lookup-tables/
      ### You can recover the "addNewData" and "readNewData" functions from GitHub:
      ### library(devtools, quietly=TRUE) #Remove the comment to load the devtools package
      ### source_gist("https://gist.github.com/dfalster/5589956") #Remove the comment to 
      ### retrieve "addNewData" and "readNewData".
    
      #### NOTE: The "expectedColumns" string in "readNewData" is modified here 
      ####       to include "source" (i.e. c("lookupVariable","lookupValue", 
      ####       "newVariable","newValue", "source")).
      ####       I included this column, because I knew I would create my lookup tables from
      ####       many sources and wanted to have a way to keep track.  This column is 
      ####       not necessary for the process to work.
  
      ####################
      ## Function descriptions for addNewData and readNewData below are
      ## modified from https://gist.github.com/dfalster/5589956
      ####################
      ##' Modifies 'data' by adding new values supplied in newDataFileName lookup table
      ##'
      ##' newDataFileName is expected to have columns 
      ##' c(lookupVariable,lookupValue,newVariable,newValue,source)
      ##' 
      ##' Within the column 'newVariable', replace values that
      ##' match 'lookupValue' within column 'lookupVariable' with the value
      ##' newValue'.  If 'lookupVariable' is NA, then replace *all* elements
      ##' of 'newVariable' with the value 'newValue'.
      ##'
      ##' Note that lookupVariable can be the same as newVariable.
      ##'
      ##' @param newDataFileName name of lookup table (e.g. title_issn_master.csv)
      ##' @param data existing data.frame (e.g. so.df)
      ##' @param allowedVars vector of permissible variable names for newVariable (e.g. ISSN)
      ##' @return modified data.frame (e.g. "title_issn_pairing.df")
      ##'

  # Define function "addNewData" (i.e. run the lines of code defining the addNewData function below.)

addNewData <- function(newDataFileName, data, allowedVars){
  
  import <- readNewData(newDataFileName, allowedVars)
  
  if( !is.null(import)){    
    for(i in seq_len(nrow(import))){  #Make replacements
      col.to <- import$newVariable[i] 
      col.from <- import$lookupVariable[i]
      if(is.na(col.from)){ # apply to whole column
        data[col.to] <- import$newValue[i]
      } else { # apply to subset
        rows <- data[[col.from]] == import$lookupValue[i]
        data[rows,col.to] <- import$newValue[i]
      }
    }   
  }      
  data
}

##' Utility function to read/process newDataFileName for addNewData
##' 
##' @param newDataFileName name of lookup table
##' @param allowedVars vector of permissible variable names for newVariable
##' @return data.frame with columns c(lookupVariable,lookupValue,newVariable,newValue,source)

### Define function "readNewData" (i.e. run the lines of code defining the readNewData function below.)

readNewData <- function(newDataFileName, allowedVars){
  
  if( file.exists(newDataFileName)){
    import <- read.csv(newDataFileName, header=TRUE, stringsAsFactors=FALSE,
                       strip.white=TRUE)
    if( nrow(import)> 0 ){
      
      #Check columns names for import are right
      expectedColumns<- c("lookupVariable","lookupValue","newVariable","newValue", "source")
      nameIsOK <-  expectedColumns %in% names(import)
      if(any(!nameIsOK))
        stop("Incorrect name in lookup table for ",
             newDataFileName, "--> ", paste(expectedColumns[!nameIsOK],
                                            collapse=", "))
      
      #Check values of newVariable are in list of allowed variables
      import$lookupVariable[import$lookupVariable == ""] <- NA
      nameIsOK <- import$newVariable %in% allowedVars
      if(any(!nameIsOK))
        stop("Incorrect name(s) in newVariable column of ",
             newDataFileName, "--> ", paste(import$newVariable[!nameIsOK],
                                            collapse=", "))
    } else {
      import <- NULL
    }
  } else {
    import <- NULL
  }
  import
}

  # DEFINE THE "ALLOWED" VARIABLE FROM THE LOOKUP TABLE
    # NOTE: This is the source of new data for the title frequency dataframe and it is these values
    # that are added to "so.df" to create the new "title_issn_pairing.df" dataframe.  In this first step
    # of adding ISSN values to journal titles, the "ALLOWED" variable is ISSN.

allowedVars <- c("ISSN") 

  # MATCH JOURNAL TITLE VARIANTS AGAINST TITLE/ISSN PAIRS IN LOOKUP TABLE
    # Use addNewData function to match journal titles in the "so.df" title frequency dataframe
    # against title/ISSN pairs in the "title_issn_master.csv" lookup table to create 
    # a new 3 column (i.e. "RefSource", "Frequency", and "ISSN") dataframe called "title_issn_pairing.df" 

title_issn_pairing.df <- addNewData("title_issn_master.csv", so.df, allowedVars)

  # VIEW "title_issn_pairing.df" AND NOTE "needsISSN" and "NA" IN THE ISSN COLUMN

View(title_issn_pairing.df)

    # NOTE: After running "addNewData" function, the ISSN column will contain 
    # one of three things: 1) the appropriate ISSN (this is what we want!), 
    # 2) the phrase "needsISSN", or 3) "NA".
    
    # NOTE: "needsISSN" is assigned when the "title_issn_master.csv" DOES HAVE A TITLE
    # that matches a "RefSource" title in the so.df dataframe, BUT DOES NOT HAVE A MATCHING ISSN.

    # This is a function of how the "title_issn_master.csv" file was created...not every 
    # title variant gets an ISSN assigned and these unassigned cells are initially empty.  
    # However, the addNewData function will not work if there are empty ISSN values.  
    # Therefore, when the "title_issn_master.csv" is created, these empty ISSN values
    # are/must be filled with something (e.g. "needsISSN"). 
    # So, when you see "needsISSN" that means THE ISSN MUST BE ADDED
    # to an existing title variant in the "title_issn_master.csv" lookup table
    # Then when the addNewData function is re-run for the so.df dataframe, that ISSN will be assigned.
    
    # NOTE: The grayed out "NA" is assigned when the so.df dataframe has a journal title variant that 
    # has no match to the RefSource title list (i.e. DOES NOT have a title 
    # (and therefore also no ISSN) that matches) in the "title_issn_master.csv" lookup 
    # table file.  This means that both THE TITLE VARIANT and THE ISSN must be 
    # added to the "title_issn_master.csv" lookup table.  Then when the addNewData 
    # function is re-run for the so.df dataframe, that ISSN will be assigned.

  ###########################################################
  ####
  ####  STEP 5C: COMBINE SOURCE TITLES - CREATE FREQUENCY DATAFRAME "issn.freq.table.df" BASED ON ISSN
  ####
  ###########################################################

  # DESCRIPTION OF THREE MAJOR PARTS:
    # 1. Convert all grayed "NA" values in the ISSN variable of the 
    #    "title_issn_pairing.df" dataframe to "Needs-title-and-ISSN".
    #     [Recall NA (i.e. empty) is generated is the ISSN column of the 
    #     "title_issn_pairing.df" dataframe for titles in so.df dataframe 
    #     that have no corresponding title variant in the in the 
    #     "title_issn_master.csv" lookup table.] 
    #
    # 2. Convert tabular data in "title_issn_pairing.df" to row data with a
    #    separate row for each "observation" based on the "Frequency" variable (e.g. a journal 
    #    with a frequency count of 10 will be represented as ten separate rows 
    #    (i.e. one for each "observation)).
    #
    # 3. Create a new frequency table/dataframe based on the number of "observations" for each ISSN.

  ### STEP 5D.1:REPLACE BLANK DATA FIELDS WITH "NA".  THEN REPLACE "NA" WITH "Needs-Title-and-ISSN"
    # At this point, empty ISSN cells in the "title_issn_pairing.df" will display 
    # with a gray "NA".
    # But these gray "NA"s just represent "empty".  For our purposes, we will actually assign
    # the value "NA" to each.
    # There is a very real possibility that your title_issn_master.csv lookup table will not have 
    # an ISSN value for every title variant in your so.df title frequency dataframe. This allows you
    # to identify those and update your title_issn_master.csv lookup table.

    # Any title without an assigned ISSN will be have the ISSN value assigned first as "NA"
    # and then all these "NA" values will be changed to "Needs-Title-and-ISSN".

title_issn_pairing.df[title_issn_pairing.df==""] <- NA

  # Replace all NA with "Needs-Title-and-ISSN" in base R

title_issn_pairing.df <- replace(title_issn_pairing.df, is.na(title_issn_pairing.df), "Needs-Title-and-ISSN")

  # Or with dplyr
  # library(dplyr)
  # title_issn_pairing.df <- title_issn_pairing.df %>% replace(., is.na(.), "Needs-Title-and-ISSN")

  ### STEP 5D.2: CONVERT TABULAR DATA IN "title_issn_pairing.df" TO ROW DATA
    # Go from the title_issn_pairing.df dataframe to a frequency table with a
    # separate row for each "observation" based on the "Frequency" variable,
    # to a dataframe listing each observation (i.e. one row per case).
  
    # An important caveat here is that your lookup table must assign an ISSN to EVERY title variant 
    # to have an absolutely accurate frequency count.  Some titles have MANY variants.  
    # So for very large issn.freq.table.df frequency tables (i.e. for an entire department or even 
    # a college), building the lookup table can take quite a while to build.  
    # The title_issn_master.csv lookup table file I created (and included here) should provide a great start.
    # It is almost a certainty that you will add to this file, so some type of file name versioning 
    # is HIGHLY RECOMMENDED.
  
    # I got this code from the Lynda.com video "R Statistics Essential Training"
    # by Barton Poulson in the section called "Converting tabular data to row data"
    # I'm using the individual steps for this transformation, starting with the
    # second step which uses a dataframe (and not the frequency table used in the example code).
      # Note that this "input" for these steps is:
      # "issn_title_pairing.df" (name of the starting dataframe)
      # "-2" (the column/variable indicting the frequency number which will be redundant information)
      # And the final output is:
      # "so.df.table5" (the name of the final dataframe containing all of the 
      # individual obervations).  This is a table with two variables (RefSource and ISSN).


so.df.table2 <- lapply(title_issn_pairing.df, function(x)rep(x, title_issn_pairing.df$Frequency)) #Use 'lapply' on title_issn_master_df dataframe and use a function which takes a variable 'x' (in this case "Freq") and repeats it the number of times is shows up in so.df.table1.  This new object (so.df.table2) is actually a list.

so.df.table3 <- as.data.frame(so.df.table2) # Converts so.df.table2 list object into a dataframe

so.df.table4 <- so.df.table3[, -2] #This creates a new dataframe where the second column/variable is removed (which is just the original frequency value and it is redundant, because we have each value in column 1 (i.e. the journal titles) repeated the appropriate number of times.)

#so.df.table4 <- so.df.table3 # If you want to keep the original title frequency data, you could just pass the so.df.table3 table into the so.df.table4 table.  Then continue to the next step.

so.df.table5 <- as.data.frame(so.df.table4) #Read in the "Large Factor" object so.df.table4 as a dataframe.

head(so.df.table5) #Look at the first several rows to confirm that the so.df.table5 dataframe seems ok.

  ###########################################################
  ####
  ####  STEP 5D: GENERATE THE "issn.freq.table.df" FREQUENCY TABLE BASED ON THE NUMBER OF OBSERVATIONS
  ###   FOR EACH ISSN IN THE "so.df.table5" DATAFRAME
  ####
  ###########################################################

    # The next step creates an ISSN frequency table (called "issn.freq.table.df") based 
    # on the number of ISSNs assigned to the title variants in the "title_issn_pairing.df" 
    # dataframe (i.e. the number of "observations" for each ISSN in the "so.df.table5" dataframe).
    # Because the ISSNs can be used as a common, unique identifier for all variants
    # of a single title, this allows you to create a more accurate frequency count of each cited source.


  # This creates a frequency table that is a dataframe.
    # Note: Using the dplyr package and the plyr package together causes problems.  
    #       If you want to use dplyr, you must load the plyr package BEFORE dplyr.

   # library(dplyr) #dplyr should have been loaded in step 0, but if not do so now.

   issn.freq.table.df <- so.df.table5 %>% plyr::count('ISSN')
   colnames(issn.freq.table.df) <- c("ISSN", "Frequency") #Renames variables for the columns
  
  ###########################################################
  ####
  ####  STEP 5E: COMBINE SOURCE TITLES - ADD STANDARDIZED JOURNAL TITLES TO "issn.freq.table.df" FREQUENCY DATAFRAME
  ####
  ###########################################################
  
  ### STEP 5E.1: ADD STANDARDIZED JOURNAL TITLES TO "issn.freq.table.df" FREQUENCY DATAFRAME
    # Now that we have a frequency table of ISSNs (issn.freq.table.df), we can use those ISSNs as
    # unique identifiers to assign a single journal title to each. This uses a second lookup table called 
    # "issn_title_master.csv" along with the addNewData and readNewData functions. 
    # This will make it easier to generate reports.  People can understand a journal
    # title more easily than an ISSN. 

    # NOTE: As with the the initial ISSN assigment in Step 
    #       Use source_gist("https://gist.github.com/dfalster/5589956") 
    #       to perform lookup table search, but this time 
    #       you will add "RefSource" data variable values.

    # NOTE: RefSource values are found in a .csv lookup table (e.g. "issn_title_master.csv", 
    #       "issn_title_freedocol.csv", etc.) and added to the issn.freq.table.df 
    #       dataframe to create the new dataframe "issn_title_pairing.df"
    
    # Note: The "issn_title_master.csv" or "issn_title_freedomcol.csv" lookup tables
    #       can't have any blank values in the lookupValue column.

  ## DEFINE THE "ALLOWED" VARIABLE (RefSource) FROM THE LOOKUP TABLE
    ### Assign "RefSource" as the allowed variable from the "issn_title_master.csv" lookup table.
  
allowedVars <- c("RefSource") #This is the variable that will be added.

  # MATCH ISSN's With ISSN/UNIQUE TITLE PAIRS IN "issn_title_master.csv" LOOKUP TABLE
    # Use addNewData function to match ISSNs in the "issn.freq.table.df" ISSN frequency dataframe
    # against ISSN/title pairs in the "issn_title_master.csv" lookup table to create 
    # a new 3 column (i.e. "RefSource", "Frequency", and "ISSN") dataframe called "title_issn_pairing.df" 
  
issn_title_pairing.df <- addNewData("issn_title_master.csv", issn.freq.table.df, allowedVars)


  ### STEP 5E.2: ORDER THE "issn_title_pairing.df" DATAFRAME BY CITATION FREQUENCY
    # Order the issn_title_pairing.df by frequency from largest to smallest
    # issn_title_pairing.df has all the data, but it is not ordered.
    # This orders issn_title_pairing.df by citation frequency in descending order

  
issn_title_pairing2.df <-issn_title_pairing.df[order(-issn_title_pairing.df$Frequency),]

  #####
  # OPTIONAL: Conditionally remove rows where $RefSource is NA.
  # In the case where we are using an ISSN to Journal Title lookup table (i.e. typically 
  # the second lookup table in this analysis) that is generated from an ISSN/Title List 
  # from a single publisher (e.g. the Elsevier Freedom Collection), there will 
  # be many ISSNs that do not have assigned Journal Titles.
  # This leaves and empty cell in the RefSource column.  These rows can be removed using na.omit.
  #####
  
  # Run only if interested in removing rows from dataframe containing "NA" values.
  issn_title_pairing2.df <- na.omit(issn_title_pairing2.df)


  
###########################################################
####
####  STEP 5F: COMBINE SOURCE TITLES - WRITE THE ORDERED DATAFRAME "issn_title_pairing2.df" TO .CSV FILE 
####
###########################################################
  
  # WRITING "Cited in" DATA TO A .CSV FILE
    # NOTE: If you use only the "write.csv" function to create the .csv file, the ISSN values are treated
    #       as numbers (and not characters).  This can be a problem is you want to use the .csv file with Excel.
    #       If this is a problem, see the optional step that follows.
   
  
write.csv(issn_title_pairing2.df, "citedIn-Example-080620.csv", 
          row.names = FALSE)

  ### [OPTIONAL] MODIFY ISSN VALUES IN DATAFRAME TO FORCE MS EXCEL TO READ ISSNs IN .CSV FILE AS CHARACTERS 
    #  NOTE: If you will use Excel to work with the CSV file, you must force the ISSN
    #  values to be read as characters (i.e. not numbers).
    #  Make sure Excel does not treat the ISSN column values as numbers.  When this happens, 
    #  if the ISSN has any 0's in the first set of the 8 postitions (e.g. 00XXXXXX) 
    #  the 0's get removed (i.e. XXXXXX).
    
  # This surrounds the ISSN values with quotes and then 
  # prepends an equal sign on to that (e.g. 00XXXXXX becomes ="00XXXXXX".)  
  # Excel reads this in as a character.

issn_title_pairing2.df$ISSN = paste0('="', issn_title_pairing2.df$ISSN, '"')
write.csv(issn_title_pairing2.df, "citedIn-Example-ISSN-as-Character-082819.csv", 
          row.names = FALSE)
