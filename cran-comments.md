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






