Clazz.declarePackage ("J.renderspecial");
Clazz.load (["J.render.ShapeRenderer", "JU.P3", "$.P3i", "$.V3"], "J.renderspecial.VectorsRenderer", ["J.shape.Shape", "JU.Vibration"], function () {
c$ = Clazz.decorateAsClass (function () {
this.pointVectorStart = null;
this.ptTemp = null;
this.pointVectorEnd = null;
this.pointArrowHead = null;
this.screenVectorStart = null;
this.screenVectorEnd = null;
this.screenArrowHead = null;
this.headOffsetVector = null;
this.diameter = 0;
this.headWidthPixels = 0;
this.vectorScale = 0;
this.vectorSymmetry = false;
this.headScale = 0;
this.drawShaft = false;
this.vibTemp = null;
this.vectorsCentered = false;
this.standardVector = true;
this.vibrationOn = false;
this.drawCap = false;
this.showModVecs = false;
Clazz.instantialize (this, arguments);
}, J.renderspecial, "VectorsRenderer", J.render.ShapeRenderer);
Clazz.prepareFields (c$, function () {
this.pointVectorStart =  new JU.P3 ();
this.ptTemp =  new JU.P3 ();
this.pointVectorEnd =  new JU.P3 ();
this.pointArrowHead =  new JU.P3 ();
this.screenVectorStart =  new JU.P3i ();
this.screenVectorEnd =  new JU.P3i ();
this.screenArrowHead =  new JU.P3i ();
this.headOffsetVector =  new JU.V3 ();
});
Clazz.overrideMethod (c$, "render", 
function () {
var vectors = this.shape;
if (!vectors.isActive) return false;
var mads = vectors.mads;
if (mads == null) return false;
var atoms = vectors.atoms;
var colixes = vectors.colixes;
var needTranslucent = false;
this.vectorScale = this.vwr.getFloat (1649410049);
if (this.vectorScale < 0) this.vectorScale = 1;
this.vectorSymmetry = this.vwr.getBoolean (603979973);
this.vectorsCentered = this.vwr.getBoolean (603979972);
this.showModVecs = this.vwr.getBoolean (603979927);
this.vibrationOn = this.vwr.tm.vibrationOn;
this.headScale = -0.2;
if (this.vectorScale < 0) this.headScale = -this.headScale;
var haveModulations = false;
for (var i = this.ms.ac; --i >= 0; ) {
var atom = atoms[i];
if (!this.isVisibleForMe (atom)) continue;
var mod = this.ms.getModulation (i);
if (this.showModVecs && !haveModulations && mod != null) haveModulations = true;
var vib = this.ms.getVibration (i, false);
if (vib == null) continue;
if (!this.transform (mads[i], atom, vib, mod)) continue;
if (!this.g3d.setC (J.shape.Shape.getColix (colixes, i, atom))) {
needTranslucent = true;
continue;
}this.renderVector (atom);
if (this.vectorSymmetry) {
if (this.vibTemp == null) this.vibTemp =  new JU.Vibration ();
this.vibTemp.setT (vib);
this.vibTemp.scale (-1);
this.transform (mads[i], atom, this.vibTemp, null);
this.renderVector (atom);
}}
if (haveModulations) for (var i = this.ms.ac; --i >= 0; ) {
var atom = atoms[i];
if (!this.isVisibleForMe (atom)) continue;
var mod = this.ms.getModulation (i);
if (mod == null) continue;
if (!this.transform (mads[i], atom, null, mod)) continue;
if (!this.g3d.setC (J.shape.Shape.getColix (colixes, i, atom))) {
needTranslucent = true;
continue;
}this.renderVector (atom);
}
return needTranslucent;
});
Clazz.defineMethod (c$, "transform", 
 function (mad, atom, vib, mod2) {
var isMod = (vib == null || vib.modDim >= 0);
var isSpin = (!isMod && vib.modDim == -2);
if (vib == null) vib = mod2;
this.drawCap = true;
if (!isMod) {
var len = vib.length ();
if (Math.abs (len * this.vectorScale) < 0.01) return false;
this.standardVector = true;
this.drawShaft = (0.1 + Math.abs (this.headScale / len) < Math.abs (this.vectorScale));
this.headOffsetVector.setT (vib);
this.headOffsetVector.scale (this.headScale / len);
}this.ptTemp.setT (atom);
var mod = atom.getModulation ();
if (this.vibrationOn && mod != null) this.vwr.tm.getVibrationPoint (mod, this.ptTemp, 1);
if (isMod) {
this.standardVector = false;
this.drawShaft = true;
mod = vib;
this.pointVectorStart.setT (this.ptTemp);
this.pointVectorEnd.setT (this.ptTemp);
if (mod.isEnabled ()) {
if (this.vibrationOn) {
this.vwr.tm.getVibrationPoint (vib, this.pointVectorEnd, NaN);
}mod.addTo (this.pointVectorStart, NaN);
} else {
mod.addTo (this.pointVectorEnd, 1);
}this.headOffsetVector.sub2 (this.pointVectorEnd, this.pointVectorStart);
var len = this.headOffsetVector.length ();
this.drawCap = (len + -0.2 > 0.001);
this.drawShaft = (len > 0.01);
this.headOffsetVector.scale (this.headScale / this.headOffsetVector.length ());
} else if (this.vectorsCentered || isSpin) {
this.standardVector = false;
this.pointVectorEnd.scaleAdd2 (0.5 * this.vectorScale, vib, this.ptTemp);
this.pointVectorStart.scaleAdd2 (-0.5 * this.vectorScale, vib, this.ptTemp);
} else {
this.pointVectorEnd.scaleAdd2 (this.vectorScale, vib, this.ptTemp);
this.screenVectorEnd.setT (this.vibrationOn ? this.tm.transformPtVib (this.pointVectorEnd, vib) : this.tm.transformPt (this.pointVectorEnd));
this.pointArrowHead.add2 (this.pointVectorEnd, this.headOffsetVector);
this.screenArrowHead.setT (this.vibrationOn ? this.tm.transformPtVib (this.pointArrowHead, vib) : this.tm.transformPt (this.pointArrowHead));
}if (!this.standardVector) {
this.screenVectorEnd.setT (this.tm.transformPt (this.pointVectorEnd));
this.screenVectorStart.setT (this.tm.transformPt (this.pointVectorStart));
if (this.drawCap) this.pointArrowHead.add2 (this.pointVectorEnd, this.headOffsetVector);
 else this.pointArrowHead.setT (this.pointVectorEnd);
this.screenArrowHead.setT (this.tm.transformPt (this.pointArrowHead));
}this.diameter = Clazz.floatToInt (mad < 0 ? -mad : mad < 1 ? 1 : this.vwr.tm.scaleToScreen (this.screenVectorEnd.z, mad));
this.headWidthPixels = this.diameter << 1;
if (this.headWidthPixels < this.diameter + 2) this.headWidthPixels = this.diameter + 2;
return true;
}, "~N,JM.Atom,JU.Vibration,J.api.JmolModulationSet");
Clazz.defineMethod (c$, "renderVector", 
 function (atom) {
if (this.drawShaft) {
if (this.standardVector) this.g3d.fillCylinderScreen (1, this.diameter, atom.sX, atom.sY, atom.sZ, this.screenArrowHead.x, this.screenArrowHead.y, this.screenArrowHead.z);
 else this.g3d.fillCylinderScreen (2, this.diameter, this.screenVectorStart.x, this.screenVectorStart.y, this.screenVectorStart.z, this.screenArrowHead.x, this.screenArrowHead.y, this.screenArrowHead.z);
}if (this.drawCap) this.g3d.fillConeScreen (2, this.headWidthPixels, this.screenArrowHead, this.screenVectorEnd, false);
}, "JM.Atom");
Clazz.defineStatics (c$,
"arrowHeadOffset", -0.2);
});
