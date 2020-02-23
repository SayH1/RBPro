headofdataframe <- function(filename, dataframe, distinct, sampletype) {


	geneDistinct <- distinct(dataframe, eval(parse(text = distinct)), .keep_all= TRUE)
	geneDistinct <- geneDistinct[,1:ncol(geneDistinct)-1]
	write.table(geneDistinct, paste("file/",filename,"_distinct.txt",sep=""), append = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)
	data <- read.delim(paste("file/",filename,"_distinct.txt",sep=""), header=T, row.names=distinct)

    # Make data matrix
    dataMatrix <- data.matrix(data, rownames.force = NA)

    # Remove non-positive values by replacing with 0
    dataMatrix[dataMatrix < 0] <- 0

    # Define group of sample using unlist function
	group <- unlist(sampletype, use.names=FALSE)

	# Design matrix (NormalvsPerturbation 0vs1)
    design <- model.matrix(~0+group)
    colnames(design) <- gsub("group", "", colnames(design))
    # print(design)

	# Construct DGEList
    x <- DGEList(counts = dataMatrix)

    # Filter genes
    keep <- filterByExpr(x, design)
    x <- x[keep,]

    # Calculate normalization factors
    x <- calcNormFactors(x)
    # x <- estimateCommonDisp(x,verbose=TRUE)

    # Assign attributes in DGEList
    # samplenames <- substring(colnames(x), 0, nchar(colnames(x)))
    # colnames(x) <- samplenames
    # group <- as.factor(unlist(sampletype, use.names=FALSE))
    # x$samples$group <- group
    # lane <- as.factor(rep(c("L004","L006","L008"), c(3,4,2)))
    # x$samples$lane <- lane
    x$samples

    # Contrast matrix (NormalvsPerturbation 1vs-1)
    contr.matrix <- makeContrasts(
       NormalvsPerturbation = Perturbation - Normal,
       levels = colnames(design))
    # contr.matrix

    # par(mfrow=c(1,2))
    v <- voom(x, plot=FALSE)
    head(v)

    vfit <- lmFit(v, design)
    vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
    efit <- eBayes(vfit)
    # plotSA(efit, main="Final model: Mean-variance trend")

    normal.vs.cancer <- topTable(efit, adjust="BH", number=nrow(dataMatrix))

    # print(head(normal.vs.cancer))
    normal.vs.cancer <- normal.vs.cancer[order(-normal.vs.cancer$t),]
    newdf <- normal.vs.cancer %>% rownames_to_column(distinct)
    write.csv(newdf,paste("file/",filename,"_Signature.txt",sep=""), sep = " ", quote=FALSE, col.names = TRUE, row.names = FALSE)

    # Quick Summary of LogCPM
    lcpm <- log(x$counts,base=2)
    # lcpm <- x$counts
    L <- mean(x$samples$lib.size) * 1e-6
    M <- median(x$samples$lib.size) * 1e-6
    c(L, M)
    summary(lcpm)

    # Examining the number of Differential expressed genes
    summary(decideTests(efit))
    edt <- decideTests(efit)
    # summary(edt)

    tfit <- treat(vfit, lfc=1)
    tdt <- decideTests(tfit)
    # summary(tdt)

    # plotMD(efit, column=1, status=edt[,1], main=colnames(efit)[1])

    # Plot interactive glimma MD plot
    glMDPlot(efit, coef=1, status=edt, main=colnames(efit)[1], counts=lcpm, groups=group, path = getwd(), folder = "static/plot", html = paste(filename, "_MD-Plot", sep=""), launch=FALSE)

    # df <- data.frame(edt)
    # ups <- which(data.frame(df) == 1, arr.ind=TRUE)
    # downs <- which(data.frame(df) == -1, arr.ind=TRUE)
    # write.table(row.names(ups), paste("file/",filename,"_Up_genes.txt",sep=""), sep = " ", quote=FALSE, row.names = FALSE, col.names = FALSE)
    # write.table(row.names(downs), paste("file/",filename,"_Down_genes.txt",sep=""), sep = " ", quote=FALSE, row.names = FALSE, col.names = FALSE)

	# Return
	return(x)
}