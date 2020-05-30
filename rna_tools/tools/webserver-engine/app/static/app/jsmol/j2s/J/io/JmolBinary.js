Clazz.declarePackage ("J.io");
Clazz.load (null, "J.io.JmolBinary", ["java.io.BufferedInputStream", "java.util.Hashtable", "JU.PT", "$.Rdr", "J.api.Interface", "JU.Logger", "JV.FileManager", "$.JmolAsyncException"], function () {
c$ = Clazz.decorateAsClass (function () {
this.fm = null;
this.pngjCache = null;
this.spardirCache = null;
this.jzu = null;
Clazz.instantialize (this, arguments);
}, J.io, "JmolBinary");
Clazz.makeConstructor (c$, 
function (fm) {
this.fm = fm;
}, "JV.FileManager");
Clazz.defineMethod (c$, "determineSurfaceTypeIs", 
function (is) {
var br;
try {
br = JU.Rdr.getBufferedReader ( new java.io.BufferedInputStream (is), "ISO-8859-1");
} catch (e) {
if (Clazz.exceptionOf (e, java.io.IOException)) {
return null;
} else {
throw e;
}
}
return this.determineSurfaceFileType (br);
}, "java.io.InputStream");
c$.getEmbeddedScript = Clazz.defineMethod (c$, "getEmbeddedScript", 
function (script) {
if (script == null) return script;
var pt = script.indexOf ("**** Jmol Embedded Script ****");
if (pt < 0) return script;
var pt1 = script.lastIndexOf ("/*", pt);
var pt2 = script.indexOf ((script.charAt (pt1 + 2) == '*' ? "*" : "") + "*/", pt);
if (pt1 >= 0 && pt2 >= pt) script = script.substring (pt + "**** Jmol Embedded Script ****".length, pt2) + "\n";
while ((pt1 = script.indexOf (" #Jmol...\u0000")) >= 0) script = script.substring (0, pt1) + script.substring (pt1 + " #Jmol...\u0000".length + 4);

if (JU.Logger.debugging) JU.Logger.debug (script);
return script;
}, "~S");
Clazz.defineMethod (c$, "getJzu", 
 function () {
return (this.jzu == null ? this.jzu = J.api.Interface.getOption ("io.JmolUtil", this.fm.vwr, "file") : this.jzu);
});
Clazz.defineMethod (c$, "getCachedPngjBytes", 
function (pathName) {
return (pathName.indexOf (".png") < 0 ? null : this.getJzu ().getCachedPngjBytes (this, pathName));
}, "~S");
Clazz.defineMethod (c$, "clearAndCachePngjFile", 
function (data) {
this.pngjCache =  new java.util.Hashtable ();
return (data == null || data[0] == null ? false : this.getJzu ().cachePngjFile (this, data));
}, "~A");
Clazz.defineMethod (c$, "spardirPut", 
function (name, bytes) {
if (this.spardirCache == null) this.spardirCache =  new java.util.Hashtable ();
this.spardirCache.put (name, bytes);
}, "~S,~A");
Clazz.defineMethod (c$, "clearPngjCache", 
function (fileName) {
if (this.pngjCache == null || fileName != null && !this.pngjCache.containsKey (fileName)) return;
this.pngjCache = null;
JU.Logger.info ("PNGJ cache cleared");
}, "~S");
Clazz.defineMethod (c$, "recachePngjBytes", 
function (fileName, bytes) {
if (this.pngjCache == null || !this.pngjCache.containsKey (fileName)) return;
this.pngjCache.put (fileName, bytes);
JU.Logger.info ("PNGJ recaching " + fileName + " (" + bytes.length + ")");
}, "~S,~A");
Clazz.defineMethod (c$, "getAtomSetCollectionOrBufferedReaderFromZip", 
function (adapter, is, fileName, zipDirectory, htParams, asBufferedReader) {
return this.getJzu ().getAtomSetCollectionOrBufferedReaderFromZip (this.fm.vwr, adapter, is, fileName, zipDirectory, htParams, 1, asBufferedReader);
}, "J.api.JmolAdapter,java.io.InputStream,~S,~A,java.util.Map,~B");
Clazz.defineMethod (c$, "spartanFileList", 
function (name, zipDirectory) {
return this.getJzu ().spartanFileList (this.fm.vwr.getJzt (), name, zipDirectory);
}, "~S,~S");
Clazz.defineMethod (c$, "determineSurfaceFileType", 
function (br) {
return this.getJzu ().determineSurfaceFileType (br);
}, "java.io.BufferedReader");
Clazz.defineMethod (c$, "getImage", 
function (vwr, fullPathNameOrBytes, echoName) {
return this.getJzu ().getImage (vwr, fullPathNameOrBytes, echoName);
}, "JV.Viewer,~O,~S");
c$.getFileReferences = Clazz.defineMethod (c$, "getFileReferences", 
function (script, fileList) {
for (var ipt = 0; ipt < JV.FileManager.scriptFilePrefixes.length; ipt++) {
var tag = JV.FileManager.scriptFilePrefixes[ipt];
var i = -1;
while ((i = script.indexOf (tag, i + 1)) >= 0) {
var s = JU.PT.getQuotedStringAt (script, i);
if (s.indexOf ("::") >= 0) s = JU.PT.split (s, "::")[1];
fileList.addLast (s);
}
}
}, "~S,JU.Lst");
c$.getManifestScriptPath = Clazz.defineMethod (c$, "getManifestScriptPath", 
function (manifest) {
if (manifest.indexOf ("$SCRIPT_PATH$") >= 0) return "";
var ch = (manifest.indexOf ('\n') >= 0 ? "\n" : "\r");
if (manifest.indexOf (".spt") >= 0) {
var s = JU.PT.split (manifest, ch);
for (var i = s.length; --i >= 0; ) if (s[i].indexOf (".spt") >= 0) return "|" + JU.PT.trim (s[i], "\r\n \t");

}return null;
}, "~S");
c$.getBufferedReaderForResource = Clazz.defineMethod (c$, "getBufferedReaderForResource", 
function (vwr, resourceClass, classPath, resourceName) {
var url;
{
}resourceName = (url == null ? vwr.vwrOptions.get ("codePath") + classPath + resourceName : url.getFile ());
if (vwr.async) {
var bytes = vwr.fm.cacheGet (resourceName, false);
if (bytes == null) throw  new JV.JmolAsyncException (resourceName);
return JU.Rdr.getBufferedReader (JU.Rdr.getBIS (bytes), null);
}return vwr.fm.getBufferedReaderOrErrorMessageFromName (resourceName, [null, null], false, true);
}, "JV.Viewer,~O,~S,~S");
Clazz.defineStatics (c$,
"JPEG_CONTINUE_STRING", " #Jmol...\0",
"PMESH_BINARY_MAGIC_NUMBER", "PM\1\0");
});
