estimate_cell_types = function(MSet){
  ## get reference dataset
  referenceRGset = get("FlowSorted.Blood.450k")
  referenceMSet = preprocessRaw(referenceRGset)
  ## generate info for combine data later
  newpd = data.frame(sampleNames = c(sampleNames(MSet), sampleNames(referenceMSet)),
                     studyIndex = c(rep("user", length(sampleNames(MSet))),
                                    rep("reference",length(sampleNames(referenceMSet)))))
  ## combine
  combine_MSet = combineArrays(MSet,referenceMSet,outType="IlluminaHumanMethylation450k")
  combine_MSet@phenoData = AnnotatedDataFrame(newpd)
  ## normalize by preprocess Quantile
  combineMSet = preprocessQuantile(combine_MSet)
  ## extract reference data and user data
  referenceMset = combineMSet[, combineMSet$studyIndex == "reference"]
  pData(referenceMset)= DataFrame(pData(referenceRGset))
  
  mSet = combineMSet[, combineMSet$studyIndex == "user"]
  pData(mSet) = as(pData(MSet),"DataFrame")
  
  compData = minfi:::pickCompProbes(referenceMset, 
                                    cellTypes = c("CD8T", "CD4T", "NK","Bcell", "Mono", "Gran"), 
                                    compositeCellType = "Blood",
                                    probeSelect = "auto")
  coefs = compData$coefEsts
  counts = minfi:::projectCellType(getBeta(mSet)[rownames(coefs), ],coefs)
  counts
}

