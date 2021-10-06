# Contributing to TETHYS

First of all, thank you for joining this project! :bow: 

## Community discussion :speech_balloon:

You can either use the _Discussions_ tab or the _issues_ to start a discussion about some problem or new feature. The rule of thumb being that issues should be relativelly self contained and problem-solving oriented while the discussion have a broader spectrum. 

## Git Flow :arrows_counterclockwise:

There are two long lived branches: `master` (as one would expect) and `TestingArea` (you can think of it as a _develop branch_). Aditionally you will find some (_'kinda'_) short-lived branches dedicated to solve or explore particular problems or physical scenarios, let's call them _feature branches_.  
While working on new aspects or features you should work on on of the existing feature branches, or if none of them adress you intended topic, create one yourself. Afterwards you can pull request/merge your addition to the `TestingArea` where seceral features are then tested and harmonized. 
Lastly those changes will be incorporated on the `master` branch. To summarize: 
<p align="center">
  Feature Branch  :arrow_right: TestingArea :arrow_right: master
</p>


## Style guide :page_with_curl:

### Semantic Versioning

Standard form of numeric *major.minor.patch* starting with the initial commit 1.0.0. Small (but relevant) bugs are considered lower level patches and new features (such as updating physical model) are minor level. Major level versions should be saved for breaking updates (like 2D implementation or parallelization)


### Internal syntax

| Type            | Style                                 | E.g.              |
| :-------------: |:-------------:                        | :-----            |
| *Macros*          | Prefix + _ + Upper-case                | MAT_PI            |
| *Functions*       | Camel case                            | InitialCondRand   |
| *Variables*       | Lower-case 3 letters code + suffix     | den_mid           |
