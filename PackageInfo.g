#############################################################################
##  
##  PackageInfo.g file for the StandardFF package.               Frank Lübeck
##  

SetPackageInfo( rec(

PackageName := "StandardFF",
Subtitle := "Standard finite fields and cyclic generators",
Version := "1.0",
Date := "13/09/2023",

Persons := [
  rec(
    FirstNames := "Frank",
    LastName := "Lübeck",
    WWWHome := "https://www.math.rwth-aachen.de/~Frank.Luebeck/",
    Email := "Frank.Luebeck@Math.RWTH-Aachen.De",
    IsAuthor := true,
    IsMaintainer := true,
    PostalAddress := "Dr.Frank Lübeck\nLehrstuhl für Algebra und Zahlentheorie\nRWTH Aachen\nPontdriesch 14/16\n52062 Aachen\nGERMANY\n",
    Place := "Aachen",
    Institution := "Lehrstuhl für Algebra und Zahlentheorie, RWTH Aachen",
  ),
],

SourceRepository := rec( Type := "git", URL := "https://github.com/frankluebeck/StandardFF" ),
IssueTrackerURL := "https://github.com/frankluebeck/StandardFF",
#SupportEmail := "",

PackageWWWHome := "https://www.math.rwth-aachen.de/~Frank.Luebeck/gap/StandardFF/",

PackageInfoURL := Concatenation( ~.PackageWWWHome, "PackageInfo.g" ),
README_URL     := Concatenation( ~.PackageWWWHome, "README.txt" ),
ArchiveURL     := Concatenation( ~.PackageWWWHome,
                                 "/", ~.PackageName, "-", ~.Version ),

ArchiveFormats := ".tar.gz",

##  Status information. Currently the following cases are recognized:
##    "accepted"      for successfully refereed packages
##    "submitted"     for packages submitted for the refereeing
##    "deposited"     for packages for which the GAP developers agreed
##                    to distribute them with the core GAP system
##    "dev"           for development versions of packages
##    "other"         for all other packages
##
Status := "deposited",
License := "GPL-3.0-or-later",

AbstractHTML   :=  "",

PackageDoc := rec(
  BookName  := "StandardFF",
  ArchiveURLSubset := ["doc"],
  HTMLStart := "doc/chap0.html",
  PDFFile   := "doc/manual.pdf",
  SixFile   := "doc/manual.six",
  LongTitle := "Standard finite fields and cyclic generators",
),

Dependencies := rec(
  GAP := ">= 4.12",
  NeededOtherPackages := [ [ "GAPDoc", ">= 1.6.5" ] ],
  SuggestedOtherPackages := [ [ "FactInt", ">= 1.6.3" ], [ "CtblLib", ">= 1.3.1" ] ],
  ExternalConditions := [ ],
),

AvailabilityTest := function()
        return true;
    end,

TestFile := "tst/testall.g",

#Keywords := [ "TODO" ],

));


