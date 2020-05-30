Clazz.declarePackage ("J.adapter.readers.xtal");
Clazz.load (["J.adapter.smarter.AtomSetCollectionReader"], "J.adapter.readers.xtal.VaspChgcarReader", ["JU.PT"], function () {
c$ = Clazz.declareType (J.adapter.readers.xtal, "VaspChgcarReader", J.adapter.smarter.AtomSetCollectionReader);
Clazz.overrideMethod (c$, "initializeReader", 
function () {
this.setSpaceGroupName ("P1");
this.setFractionalCoordinates (true);
});
Clazz.overrideMethod (c$, "checkLine", 
function () {
var atomSym = this.getTokens ();
var scale = this.parseFloatStr (this.rd ());
var unitCellData =  Clazz.newFloatArray (9, 0);
this.fillFloatArray (null, 0, unitCellData);
for (var i = 0; i < 9; i++) unitCellData[i] *= scale;

this.addPrimitiveLatticeVector (0, unitCellData, 0);
this.addPrimitiveLatticeVector (1, unitCellData, 3);
this.addPrimitiveLatticeVector (2, unitCellData, 6);
var tokens = JU.PT.getTokens (this.rd ());
var atomCounts =  Clazz.newIntArray (tokens.length, 0);
for (var i = tokens.length; --i >= 0; ) atomCounts[i] = this.parseIntStr (tokens[i]);

if (atomSym.length != atomCounts.length) atomSym = null;
this.rd ();
for (var i = 0; i < atomCounts.length; i++) for (var j = atomCounts[i]; --j >= 0; ) this.addAtomXYZSymName (JU.PT.getTokens (this.rd ()), 0, (atomSym == null ? "Xx" : atomSym[i]), null);


this.continuing = false;
return false;
});
Clazz.overrideMethod (c$, "finalizeSubclassReader", 
function () {
this.applySymmetryAndSetTrajectory ();
});
});
