# v1.2.0

* This is a new release.

## Initial Submission

### Test environments
* local R installation, R 4.0.2
* ubuntu 16.04 (on travis-ci), R 4.0.2
* win-builder (devel)

### R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Jui-Hua Hsieh <juihua.hsieh@gmail.com>'

New submission

Found the following (possibly) invalid URLs:
  URL: https://doi.org/10.1093/toxsci/kfy258
    From: man/Rcurvep.Rd
    Status: Error
    Message: libcurl error code 56:
      	Send failure: Connection was reset
      	
The link is valid.

## Resubmission

Response to the review by Gregor Seyer

#--------

REQUEST:
> Please reduce the length of the title to less than 65 characters.

RESPONSE:
The title is now shortened. 

#---------

REQUEST: 
> If there are references describing the methods in your package, please
add these in the description field of your DESCRIPTION file in the form
authors (year) <doi:...>
authors (year) <arXiv:...>
authors (year, ISBN:...)
or if those are not available: <https:...>
with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for
auto-linking.
(If you want to add a title as well please put it in quotes: "Title")

RESPONSE:
The DESCRIPTION file is updated with references.

#---------

REQUEST:
>\dontrun{} should only be used if the example really cannot be executed
(e.g. because of missing additional software, missing API keys, ...) by
the user. That's why wrapping examples in \dontrun{} adds the comment
("# Not run:") as a warning for the user.
Does not seem necessary.
Please unwrap the examples if they are executable in < 5 sec, or replace
\dontrun{} with \donttest{}.

RESPONSE:
All the dontrun blocks are replaced by the donttest blocks.

#---------


# v1.2.1

## check results

### Test environments

#### local R installation, R 4.2.1

0 errors, 0 warnings, 0 notes

### Rhub  Windows Server 2022, R-devel, 64 bit

0 errors, 0 warnings, 1 note


❯ checking for detritus in the temp directory ... NOTE
  Found the following files/directories:
    'lastMiKTeXException'
    
Response:it seems to a bug/crash in miktex, suggested to be ignored ([R-hub issue #503](https://github.com/r-hub/rhub/issues/503)

### Rhub  Ubuntu Linux 20.04.1 LTS, R-release, GCC

0 errors, 0 warnings, 1 note

* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Jui-Hua Hsieh <juihua.hsieh@gmail.com>’

Found the following (possibly) invalid DOIs:
  DOI: 10.1093/toxsci/kfy258
    From: DESCRIPTION
    Status: Forbidden
    Message: 403
    
Response:the DOI link is working but it not open access.  
    
### Rhub Fedora Linux, R-devel, clang, gfortran

0 errors, 0 warnings, 1 note

* checking HTML version of manual ... NOTE
Skipping checking HTML validation: no command 'tidy' found

Response: it seems to the environment variable setting in R and the tidy command issue ([R release news](https://cran.r-project.org/doc/manuals/r-release/NEWS.html)


### winbuilder

0 errors, 0 warnings, 0 notes





