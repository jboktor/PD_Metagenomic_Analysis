#!/bin/bash
# adapted from: https://github.com/jladau/Nestedness.git

helpFunction()
{
   echo ""
   echo "Usage: $0 -a parameterA -b parameterB -c parameterC -d parameterD"
   echo -e "\t-a Description of what is parameterA"
   echo -e "\t-b Description of what is parameterB"
   echo -e "\t-c Description of what is parameterC"
   echo -e "\t-d Description of what is parameterD"
   exit 1 # Exit script after printing help
}

while getopts "a:b:c:d:" opt
do
   case "$opt" in
      a ) parameterA="$OPTARG" ;;
      b ) parameterB="$OPTARG" ;;
      c ) parameterC="$OPTARG" ;;
      d ) parameterD="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$parameterA" ] || [ -z "$parameterB" ] || [ -z "$parameterC" ]|| [ -z "$parameterD" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

# Begin script in case all parameters are correct
echo "$parameterA"
echo "$parameterB"
echo "$parameterC"
echo "$parameterD"

sIODir=$parameterA
sOutputDir=$sIODir$parameterB
sBiomPath=$sIODir$parameterC

# sIODir=$HOME/Documents/Research/Java/Distribution/Nestedness
# sOutputDir=$sIODir/test-output
# sBiomPath=$sIODir/test-data/test-data.biom
sJavaDir=$sIODir/bin
sTaxonRank=otu
sNullModel=equiprobablefixed
sAxis=sample
sComparisonMode=betweeneachpairoftypes

mkdir -p $sOutputDir
cd $sOutputDir

#loading sample subset
#java -cp $sJavaDir/Utilities.jar edu.ucsf.BIOM.PrintMetadata.PrintMetadataLauncher --help > $sIODir/doc/Utilities.edu.ucsf.BIOM.PrintMetadata.PrintMetadataLauncher.txt
java -cp $sJavaDir/Utilities.jar edu.ucsf.BIOM.PrintMetadata.PrintMetadataLauncher \
	--sBIOMPath=$sBiomPath \
	--sOutputPath=$sOutputDir/temp.0.csv \
	--sAxis=sample
cut -d\, -f1 temp.0.csv | head -n 5000 > samples-to-keep.csv

#making nestedness graph
#java -cp $sJavaDir/Autocorrelation.jar edu.ucsf.Nestedness.Grapher.GrapherLauncher --help > $sIODir/doc/Autocorrelation.edu.ucsf.Nestedness.Grapher.GrapherLauncher.txt
java -cp $sJavaDir/Autocorrelation.jar edu.ucsf.Nestedness.Grapher.GrapherLauncher \
	--sSamplesToKeepPath=$sOutputDir/samples-to-keep.csv \
	--sBIOMPath=$sBiomPath \
	--bNormalize=false \
	--sTaxonRank=$sTaxonRank \
	--sOutputPath=$sOutputDir/graphs.csv \
	--rgsSampleMetadataFields=latitude

#loading comparisons
#java -cp $sJavaDir/Autocorrelation.jar edu.ucsf.Nestedness.ComparisonSelector.ComparisonSelectorLauncher --help > $sIODir/doc/Autocorrelation.edu.ucsf.Nestedness.ComparisonSelector.ComparisonSelectorLauncher.txt
java -Xmx5g -cp $sJavaDir/Autocorrelation.jar edu.ucsf.Nestedness.ComparisonSelector.ComparisonSelectorLauncher \
	--sBIOMPath=$sBiomPath \
	--sOutputPath=$sOutputDir/comparisons.csv \
	--bNormalize=false \
	--sTaxonRank=$sTaxonRank \
	--sMetadataField=$parameterD \
	--iRandomSeed=1234 \
	--sComparisonMode=$sComparisonMode \
	--iNestednessPairs=250 \
	--sSamplesToKeepPath=$sOutputDir/samples-to-keep.csv \
	--sNestednessAxis=$sAxis \
	--iPrevalenceMinimum=1

#running statistics
#java -cp $sJavaDir/Autocorrelation.jar edu.ucsf.Nestedness.Calculator.CalculatorLauncher --help > $sIODir/doc/Autocorrelation.edu.ucsf.Nestedness.Calculator.CalculatorLauncher.txt
java -cp $sJavaDir/Autocorrelation.jar edu.ucsf.Nestedness.Calculator.CalculatorLauncher \
	--sBIOMPath=$sBiomPath \
	--sOutputPath=$sOutputDir/statistics.csv \
	--bNormalize=false \
	--sTaxonRank=$sTaxonRank \
	--sComparisonsPath=$sOutputDir/comparisons.csv \
	--iNullModelIterations=10000 \
	--bOrderedNODF=false \
	--sNestednessAxis=$sAxis \
	--sSamplesToKeepPath=$sOutputDir/samples-to-keep.csv \
	--sNestednessNullModel=$sNullModel \
	--iPrevalenceMinimum=1 \
	--bSimulate=false

# # #cleaning up
# # rm $sOutputDir/temp.*.*



