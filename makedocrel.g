##  this creates the documentation, needs: GAPDoc package, latex, pdflatex,
##  mkindex

# substitute this by path to main GAP directory, if this package is not
# in standard location
if IsBound(pathtoroot) then
  relpath := pathtoroot;
else
  relpath:="../../..";
fi;

hasext := function(nam, exts)
  local l;
  l := Length(nam);
  return ForAny(exts, ex-> l >= Length(ex) and 
                           nam{[l-Length(ex)+1..l]} = ex); 
end;
cont := DirectoryContents("lib");
cont := Filtered(cont, nam-> hasext(nam, [".g", ".gd", ".gi"]));
files := List(cont, nam-> Concatenation("../lib/", nam));
LoadPackage("GAPDoc");
tree := MakeGAPDocDoc("doc", "StandardFF", files, "StandardFF", 
                                           relpath, "MathJax");
GAPDocManualLab("StandardFF");
CopyHTMLStyleFiles("doc");

# further checks
#   add after release of GAPDoc 1.6.4
##  Print("Validating XML ...\n");
##  v := ValidateGAPDoc([tree.input, tree.inputorigins]);
Print("Checking all manual examples ...\n");
ex:=ExtractExamplesXMLTree(tree,"Chapter");
RunExamples(ex);
#QUIT;
