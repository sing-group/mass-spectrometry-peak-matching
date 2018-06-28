# MIT License
# 
# Copyright (c) 2018 SING Group

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# 

# Peak representation
# A peak is represented by its mass, intensity and an id
peak <- setClass("peak", 
	representation(
		mass="numeric",
		int="numeric",
		id="character"
));

setMethod("show", "peak",
    function(object){
		cat(paste(" m/z: ", getMass(object), "\n", ""));
		cat(paste("intensity: ", getIntensity(object), "\n", ""));
	}
)

getId <- function(peak){
	attributes(peak)$id
}

getMass <- function(peak){
	attributes(peak)$mass
}

getIntensity <- function(peak){
	attributes(peak)$int
}

massesToPeakList <- function(masses, id){
	toret <- vector();
	for(i in 1:length(masses)){
		toret <- c(toret, peak(mass=masses[i],int=1,id=id))
	}
	toret
}

# Merge two ordered lists of peak objects generating a new ordered list 
mergePeakLists <- function(x,y){
	toret <- vector();
	xI <- 1;
	yI <- 1;
	while (xI <= length(x)
		&& yI <= length(y))
	{
		if(getMass(x[[xI]]) < getMass(y[[yI]])){
			toret <- c(toret, x[[xI]]);
			xI <- xI+1;
		} else {
			toret <- c(toret, y[[yI]]);
			yI <- yI+1;
		}
	}

	if (xI > length(x)){
		toret <- c(toret,y[yI:length(y)]);
	} else if (yI > length(y)){
		toret <- c(toret,x[xI:length(x)]);
	}
}


match <- function(type, tolerance, referenceMass, mass){
	dif <- abs(referenceMass - mass);
	toret <- FALSE;

	if (type == "absolute"){
		if(dif <= tolerance){
			toret <- TRUE
		} else { toret <- FALSE }
	}
	if (type == "relative"){
		tolerance <- referenceMass * tolerance;
		if(dif <= tolerance){
			toret <- TRUE
		} else { toret <- FALSE }
	}
	if (type=="ppm"){	
		tolerance <- (referenceMass * tolerance)  / 1000000 ;
		if(dif <= tolerance){
			toret <- TRUE
		} else { toret <- FALSE }
	}
	toret
}

# Peak alignment algorithm
unifyMasses <- function(peakList, toleranceType, tolerance, verbose){
	alignedPeakList <- list();
	while(length(peakList) > 0){
		firstPeak <- peakList[[1]];
		centroid <- attributes(firstPeak)$mass;
		peakList[[1]] <- NULL;

		group <- list();
		addedSpectra <- vector();

		group <- c(group, firstPeak);
		addedSpectra <- c(addedSpectra, attributes(firstPeak)$id);
		
		finished <- FALSE;
		while(!finished
			&& length(peakList) > 0){
			
			currentPeak <- peakList[[1]];
			firstMatch <- match(toleranceType, tolerance, centroid, attributes(firstPeak)$mass);
			currentMatch <- match(toleranceType, tolerance, centroid, attributes(currentPeak)$mass);
			if (firstMatch && currentMatch){
				peakList[[1]] <- NULL;
				if (attributes(currentPeak)$id %in% addedSpectra){
					warning(paste("Warn. Ignoring peak ", getMass(currentPeak), " of spectrum ", getId(currentPeak),""));
				} else {
					group <- c(group, currentPeak);
					addedSpectra <- c(addedSpectra, attributes(currentPeak)$id);
					avgMass <- 0;
					for(i in 1:length(group)){		
						avgMass <- avgMass + attributes(group[[i]])$mass;
					}
					avgMass <- avgMass / length(group);
					centroid <- avgMass;
				}
			} else {
				finished <- TRUE;
			}
		}
		if (verbose == TRUE){
			for(i in 1:length(group)){		
				alignedPeakList <- c(alignedPeakList, peak(mass=centroid,int=getIntensity(group[[i]]),id=getId(group[[i]])));
			}
		} else {
			for(i in 1:length(group)){		
				alignedPeakList <- c(alignedPeakList, peak(mass=centroid,int=getIntensity(group[[i]]),id=getId(group[[i]])));
			}
		}
	}
	alignedPeakList
}

# Main function to use the forward algorithm
binPeaks.forward <- function(spectra, toleranceType, tolerance, verbose){

	# Default parameters
	if (missing(toleranceType)) toleranceType <- "ppm";
	if (missing(tolerance)) tolerance <- 250;
	if (missing(verbose)) verbose <- FALSE;
	
	# Add peaks of each spectrum to an ordered list
	peakList <- vector();
	for(i in 1:length(spectra)){
		peakList <- mergePeakLists(peakList, massesToPeakList(spectra[[i]], as.character(i)));
	}
		
	# Align them using unifyMasses
	alignedPeakList <- unifyMasses(peakList, toleranceType, tolerance, verbose);
	
	# Reconstruct the spectra with the aligned peaks
	alignedSpectra <- rep(list(vector()), length(spectra));
	for(i in 1:length(alignedPeakList)){
		currentPeak <- alignedPeakList[[i]];
		alignedSpectra[[as.numeric(getId(currentPeak))]] <- c(alignedSpectra[[as.numeric(getId(currentPeak))]], getMass(currentPeak));
	}
	alignedSpectra
}

