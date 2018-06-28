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

source("alignment.R")

s1 <- read.csv("data/A-HCCA.1-1.csv")
s2 <- read.csv("data/A-HCCA.1-2.csv")
s3 <- read.csv("data/A-HCCA.1-3.csv")
s4 <- read.csv("data/A-HCCA.1-4.csv")
s5 <- read.csv("data/A-HCCA.1-5.csv")
samples <- list(s1[,1], s2[,1], s3[,1], s4[,1], s5[,1])

alignedSamples <- binPeaks.forward(samples, "ppm", 250, verbose = TRUE)

print(alignedSamples)