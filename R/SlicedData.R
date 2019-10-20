SlicedData <- setRefClass( "SlicedData",
                           fields = list(
                             dataFrame = "data.frame",
                             dataEnv = "environment",
                             nSlices1 = "numeric",
                             rowNameSlices = "list",
                             columnNames = "character",
                             fileDelimiter = "character",
                             fileSkipColumns = "numeric",
                             fileSkipRows = "numeric",
                             fileSliceSize = "numeric",
                             fileOmitCharacters = "character"
                           ),
                           methods = list(
                             initialize = function( mat = NULL ) {
                               dataEnv <<- new.env(hash = TRUE, size = 29L);
                               nSlices1 <<- 0L;
                               if(!is.null(mat)) {
                                 CreateFromMatrix(mat);
                               }
                               fileSliceSize <<- 1000;
                               fileDelimiter <<- "\t";
                               fileSkipColumns <<- 1L;
                               fileSkipRows <<- 1L;
                               fileOmitCharacters <<- "NA"
                               return(invisible(.self));
                             },
                             CreateFromMatrix = function( mat ) {
                               stopifnot( class(mat) == "matrix" );
                               setSliceRaw( 1L ,mat );
                               rns = rownames( mat, do.NULL = FALSE);
                               #if( is.null(rns) ) {
                               #	rns = paste( "Row_",(1:nrow(mat)), sep="" );
                               #}
                               rowNameSlices <<- list(rns);
                               cns = colnames( mat, do.NULL = FALSE );
                               #if( is.null(cns) ){
                               #	cns = paste( "Col_",(1:ncol(mat)), sep="" );
                               #}
                               columnNames <<- cns;
                               return(invisible(.self));
                             },
                             getSlice = function(sl) {
                               value = get(paste(sl), dataEnv);
                               if( is.raw(value) ) {
                                 storage.mode(value) = "double";
                                 value[value == 255] = NA;
                               }
                               return( value  )
                             },
                             getSliceRaw = function(sl) {
                               return( get(paste(sl), dataEnv) )
                             },
                             setSliceRaw = function(sl, value) {
                               assign( paste(sl), value, dataEnv )
                               if( nSlices1 < sl ) {
                                 nSlices1 <<- sl;
                               }
                             },
                             setSlice = function(sl, value) {
                               if( length(value) > 0 ) {
                                 if( all(as.integer(value) == value, na.rm = TRUE) ) {
                                   if( (min(value, na.rm = TRUE) >= 0 ) &&
                                       (max(value, na.rm = TRUE) < 255) )
                                   {
                                     nv = value;
                                     suppressWarnings({storage.mode(nv) = "raw"});
                                     nv[ is.na(value)] = as.raw(255);
                                     value = nv;
                                   } else {
                                     storage.mode(value) = "integer";
                                   }
                                 }
                               }
                               setSliceRaw(sl, value);
                             },
                             nSlices = function() {
                               return( nSlices1 );
                             },
                             LoadFile = function(filename, skipRows = NULL, skipColumns = NULL, sliceSize = NULL, omitCharacters = NULL, delimiter = NULL, rowNamesColumn = 1) {
                               if( !is.null(skipRows) ) {
                                 fileSkipRows <<- skipRows;
                               }
                               if( !is.null(skipColumns) ) {
                                 fileSkipColumns <<- skipColumns;
                               }
                               if( !is.null(omitCharacters) ) {
                                 fileOmitCharacters <<- omitCharacters;
                               }
                               if( !is.null(sliceSize) ) {
                                 fileSliceSize <<- sliceSize;
                               }
                               if( !is.null(delimiter) ) {
                                 fileDelimiter <<- delimiter;
                               }
                               stopifnot( (fileSkipColumns == 0) || (rowNamesColumn <= fileSkipColumns) )
                               stopifnot( (fileSkipColumns == 0) || (rowNamesColumn >= 1) )

                               fid = file(description = filename, open = "rt", blocking = FALSE, raw = FALSE)
                               # clean object if file is open
                               Clear();
                               lines = readLines(con = fid, n = max(fileSkipRows,1L), ok = TRUE, warn = TRUE)
                               line1 = tail(lines,1);
                               splt = strsplit(line1, split = fileDelimiter, fixed = TRUE);
                               if( fileSkipRows > 0L ) {
                                 columnNames <<- splt[[1]]; # [ -(1:fileSkipColumns) ];
                               } else {
                                 seek(fid, 0)
                               }

                               rm( lines, line1, splt );

                               rowNameSlices <<- vector("list", 15);

                               curSliceId = 0L;
                               repeat
                               {
                                 # preallocate names and data
                                 if(length(rowNameSlices) < curSliceId) {
                                   rowNameSlices[[2L*curSliceId]] <<- NULL;
                                 }
                                 curSliceId = curSliceId + 1L;

                                 # read sliceSize rows
                                 rowtag = vector("character",fileSliceSize);
                                 rowvals = vector("list",fileSliceSize);
                                 for(i in 1:fileSliceSize) {
                                   temp = "";
                                   if( fileSkipColumns > 0L ) {
                                     temp = scan(file = fid, what = character(), n = fileSkipColumns, quiet = TRUE,sep = fileDelimiter);
                                   }
                                   rowtag[i] = temp[rowNamesColumn];#paste(temp,collapse=" ");
                                   rowvals[[i]] = scan(file = fid, what = double(), nlines = 1, quiet = TRUE, sep = fileDelimiter, na.strings = fileOmitCharacters);
                                   if( length(rowvals[[i]]) == 0L ) {
                                     if(i==1L) {
                                       rowtag = matrix(0, 0, 0);
                                       rowvals = character(0);
                                     } else 	{
                                       rowtag  = rowtag[  1:(i-1) ];
                                       rowvals = rowvals[ 1:(i-1) ];
                                     }
                                     break;
                                   }
                                 }
                                 if( length(rowtag) == 0L ) {
                                   curSliceId = curSliceId - 1L;
                                   break;
                                 }
                                 rowNameSlices[[curSliceId]] <<- rowtag;
                                 data = c(rowvals, recursive = TRUE);
                                 dim(data) = c(length(rowvals[[1]]), length(rowvals));
                                 data = t(data);
                                 setSlice(curSliceId, data);
                                 if( length(rowtag) < fileSliceSize ) {
                                   break;
                                 }
                                 numtxt = formatC(curSliceId*fileSliceSize, big.mark=",", format = "f", digits = 0)
                                 #cat( "Rows read: ", numtxt, "\n");
                                 flush.console()
                               }
                               close(fid)
                               if( fileSkipRows == 0 ) {
                                 columnNames <<- paste("Col_", (1:nCols()), sep="");
                               } else {
                                 columnNames <<- tail(columnNames, ncol(getSliceRaw(1)));
                               }
                               if( fileSkipColumns == 0 ) {
                                 cnt = 0L;
                                 for( sl in 1:nSlices() ) {
                                   nr = length(getSliceRaw(sl));
                                   rowNameSlices[[sl]] <<- paste("Row_",cnt + (1:nr),sep="");
                                   cnt = cnt + nr;
                                 }
                               }
                               rowNameSlices <<- rowNameSlices[1:curSliceId];
                               #cat("Rows read: ", nRows(), " done.\n");
                               return(invisible(.self));
                             },
                             SaveFile = function(filename) {
                               if( nSlices() == 0 ) {
                                 cat("No data to save");
                                 return();
                               }
                               fid = file(filename,"wt");
                               for( sl in 1:nSlices() ) {
                                 z = getSlice(sl);
                                 rownames(z) = rowNameSlices[[sl]];
                                 colnames(z) = columnNames;
                                 write.table(z, file = fid, sep = "\t",
                                             col.names = (if(sl == 1){NA}else{FALSE}));
                               }
                               close(fid);
                             },
                             nRows = function() {
                               s = 0L;
                               for(sl in .seq(1,nSlices())) {
                                 s = s + nrow(getSliceRaw(sl));
                               }
                               return( s )
                             },
                             nCols = function() {
                               if( nSlices() == 0L ) {
                                 return(0L);
                               } else {
                                 return( ncol(getSliceRaw(1L)) )
                               }
                             },
                             Clear = function() {
                               for( sl in .seq(1,nSlices()) ) {
                                 rm(list = paste(sl), envir = dataEnv)
                               }
                               nSlices1 <<- 0L;
                               rowNameSlices <<- list();
                               columnNames <<- character();
                               return(invisible(.self));
                             },
                             IsCombined = function() {
                               return( nSlices() <= 1L );
                             },
                             GetAllRowNames = function() {
                               return( c(rowNameSlices, recursive=TRUE) );
                             },
                             GetNRowsInSlice = function(sl) {
                               return( length( rowNameSlices[[sl]] ) );
                             },
                             SetNanRowMean = function() {
                               if( (nCols() == 0L) ) {
                                 return(invisible(.self));
                               }
                               for( sl in .seq(1,nSlices()) ) {
                                 slice = getSlice(sl);
                                 if( any(is.na(slice)) ) {
                                   rowmean = rowMeans(slice, na.rm = TRUE);
                                   rowmean[is.na(rowmean)] = 0L;
                                   for( j in which(!complete.cases(slice)) ) {
                                     where1 = is.na(slice[j, ]);
                                     slice[j, where1] = rowmean[j];
                                   }
                                   setSlice(sl, slice);
                                 }
                               }
                               return(invisible(.self));
                             },
                             RowStandardizeCentered = function() {
                               for(sl in .seq(1,nSlices()) ) {
                                 slice = getSlice(sl);
                                 div = sqrt( rowSums(slice^2) );
                                 div[ div == 0 ] = 1;
                                 setSlice(sl, slice/div);
                               }
                               return(invisible(.self));
                             },
                             CombineInOneSlice = function() {
                               if( nSlices() <= 1L ) {
                                 return(invisible(.self));
                               }
                               nc = nCols();
                               nr = nRows();
                               datatypes = c("raw","integer","double");
                               datafuns = c(as.raw, as.integer, as.double);
                               datatype = character(nSlices());
                               for(sl in 1:nSlices()) {
                                 datatype[sl] = typeof(getSliceRaw(sl));
                               }
                               mch = max(match(datatype,datatypes,nomatch = length(datatypes)));
                               datafun = datafuns[[mch]];
                               newData = matrix(datafun(0), nrow = nr, ncol = nc);
                               offset = 0;
                               for(sl in 1:nSlices()) {
                                 if(mch==1) {
                                   slice = getSliceRaw(sl);
                                 } else {
                                   slice = getSlice(sl);
                                 }
                                 newData[ offset + (1:nrow(slice)),] = datafun(slice);
                                 setSlice(sl, numeric());
                                 offset = offset + nrow(slice);
                               }

                               nSlices1 <<- 1L;
                               setSliceRaw(1L, newData);
                               rm(newData);

                               newrowNameSlices = GetAllRowNames();
                               rowNameSlices <<- list(newrowNameSlices)
                               return(invisible(.self));
                             },
                             ResliceCombined = function(sliceSize = -1) {
                               if( sliceSize > 0L ) {
                                 fileSliceSize <<- sliceSize;
                               }
                               if( fileSliceSize <= 0 ) {
                                 fileSliceSize <<- 1000;
                               }
                               if( IsCombined() ) {
                                 nRows1 = nRows();
                                 if(nRows1 == 0L) {
                                   return(invisible(.self));
                                 }
                                 newNSlices = floor( (nRows1 + fileSliceSize - 1)/fileSliceSize );
                                 oldData = getSliceRaw(1L);
                                 #oldNames = rowNameSlices[[1]];
                                 newNameslices = vector("list",newNSlices)
                                 for( sl in 1:newNSlices ) {
                                   range = (1+(sl-1)*fileSliceSize) : (min(nRows1,sl*fileSliceSize));
                                   newpart = oldData[range, ,drop = FALSE];
                                   if( is.raw(oldData) ) {
                                     setSliceRaw( sl, newpart);
                                   } else {
                                     setSlice( sl, newpart);
                                   }
                                   newNameslices[[sl]] = rowNameSlices[[1]][range];
                                 }
                                 rowNameSlices <<- newNameslices ;
                               } else {
                                 stop("Reslice of a sliced matrix is not supported yet. Use CombineInOneSlice first.");
                               }
                               return(invisible(.self));
                             },
                             Clone = function() {
                               clone = SlicedData$new();
                               for(sl in .seq(1,nSlices()) ) {
                                 clone$setSliceRaw(sl,getSliceRaw(sl));
                               }
                               clone$rowNameSlices = rowNameSlices;
                               clone$columnNames = columnNames;
                               clone$fileDelimiter = fileDelimiter;
                               clone$fileSkipColumns = fileSkipColumns;
                               clone$fileSkipRows = fileSkipRows;
                               clone$fileSliceSize = fileSliceSize;
                               clone$fileOmitCharacters = fileOmitCharacters;
                               return( clone );
                             },
                             RowMatrixMultiply = function(multiplier) {
                               for(sl in .seq(1,nSlices()) ) {
                                 setSlice(sl, getSlice(sl) %*% multiplier);
                               }
                               return(invisible(.self));
                             },
                             ColumnSubsample = function(subset) {
                               for(sl in .seq(1,nSlices()) ) {
                                 setSliceRaw(sl, getSliceRaw(sl)[ ,subset, drop = FALSE]);
                               }
                               columnNames <<- columnNames[subset];
                               return(invisible(.self));
                             },
                             RowReorderSimple = function(ordr) {
                               # had to use an inefficient and dirty method
                               # due to horrible memory management in R
                               if( (typeof(ordr) == "logical") && all(ordr) ) {
                                 return(invisible(.self));
                               }
                               if( (length(ordr) == nRows()) && all(ordr == (1:length(ordr))) ) {
                                 return(invisible(.self));
                               }
                               CombineInOneSlice();
                               gc();
                               setSliceRaw( 1L, getSliceRaw(1L)[ordr, ] );
                               rowNameSlices[[1]] <<- rowNameSlices[[1]][ordr];
                               gc();
                               ResliceCombined();
                               gc();
                               return(invisible(.self));
                             },
                             RowReorder = function(ordr) {
                               # transform logical into indices
                               if( typeof(ordr) == "logical" ) {
                                 if( length(ordr) == nRows() ) {
                                   ordr = which(ordr);
                                 } else {
                                   stop("Parameter \"ordr\" has wrong length")
                                 }
                               }
                               ## first, check that anything has to be done at all
                               if( (length(ordr) == nRows()) && all(ordr == (1:length(ordr))) ) {
                                 return(invisible(.self));
                               }
                               ## check bounds
                               #if( (min(ordr) < 1) || (max(ordr) > nRows()) ) {
                               #	stop("Parameter \"ordr\" is out of bounds");
                               #}
                               ## slice the data into individual rows
                               all_rows = vector("list", nSlices())
                               for( i in 1:nSlices() ) {
                                 slice = getSliceRaw(i)
                                 all_rows[[i]] = split(slice, 1:nrow(slice))
                                 setSliceRaw(i,numeric())
                               }
                               gc();
                               all_rows = unlist(all_rows, recursive=FALSE, use.names = FALSE);
                               ## Reorder the rows
                               all_rows = all_rows[ordr];
                               ## get row names
                               all_names = GetAllRowNames();
                               ## erase the set
                               rowNameSlices <<- list();
                               ## sort names
                               all_names = all_names[ordr];
                               ##
                               ## Make slices back
                               nrows = length(all_rows);
                               nSlices1 <<- as.integer((nrows+fileSliceSize-1)/fileSliceSize);
                               ##cat(nrows, " ", nSlices1);
                               rowNameSlices1 = vector("list", nSlices1);
                               for( i in 1:nSlices1 ) {
                                 fr = 1 + fileSliceSize*(i-1);
                                 to = min( fileSliceSize*i, nrows);

                                 subset = all_rows[fr:to];
                                 types = unlist(lapply(subset,typeof));
                                 israw = (types == "raw")
                                 if(!all(israw == israw[1])) {
                                   # some raw and some are not
                                   subset = lapply(subset, function(x){if(is.raw(x)){x=as.integer(x);x[x==255] = NA;return(x)}else{return(x)}});
                                 }
                                 subset = unlist(subset);
                                 dim(subset) = c( length(all_rows[[fr]]) , to - fr + 1)
                                 #subset = matrix(subset, ncol = (to-fr+1));
                                 if(is.raw(subset)) {
                                   setSliceRaw(i, t(subset));
                                 } else {
                                   setSlice(i, t(subset));
                                 }
                                 rowNameSlices1[[i]] = all_names[fr:to];
                                 all_rows[fr:to] = 0;
                                 all_names[fr:to] = 0;
                               }
                               rowNameSlices <<- rowNameSlices1;
                               gc();
                               return(invisible(.self));
                             },
                             RowRemoveZeroEps = function(){
                               for(sl in .seq(1,nSlices()) ) {
                                 slice = getSlice(sl);
                                 amean = rowMeans(abs(slice));
                                 remove = (amean < .Machine$double.eps*nCols());
                                 if(any(remove)) {
                                   rowNameSlices[[sl]] <<- rowNameSlices[[sl]][!remove];
                                   setSlice(sl, slice[!remove, , drop = FALSE]);
                                 }
                               }
                               return(invisible(.self));
                             },
                             FindRow = function(rowname) {
                               for(sl in .seq(1,nSlices()) ) {
                                 mch = match(rowname,rowNameSlices[[sl]], nomatch = 0);
                                 if( mch > 0 )
                                 {
                                   row = getSlice(sl)[mch[1], , drop=FALSE];
                                   rownames(row) = rowname;
                                   colnames(row) = columnNames;
                                   return( list(slice = sl, item = mch, row = row) );
                                 }
                               }
                               return( NULL );
                             },
                             show = function() {
                               #cat("SlicedData object. For more information type: ?SlicedData\n");
                               #cat("Number of columns:", nCols(), "\n");
                               #cat("Number of rows:", nRows(), "\n");
                               #cat("Data is stored in", nSlices(), "slices\n");
                               if(nCols()>0) {
                                 z = getSlice(1L);
                                 if(nrow(z)>0) {
                                   z = z[1:min(nrow(z),10L), 1:min(ncol(z),10L), drop = FALSE];
                                   rownames(z) = rowNameSlices[[1]][1:nrow(z)];
                                   colnames(z) = columnNames[1:ncol(z)];
                                   #cat("Top left corner of the first slice (up to 10x10):\n");
                                   methods:::show(z)
                                 }
                               }
                             }
                           ))





