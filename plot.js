////////////////////////////////////////////////////////////////////////////////////////////////////////
////// created by Jan £ukasz Górski, janlgorski@gmail.com, http://mathmed.blox.pl 100% opensource //////
////////////////////////////////////////////////////////////////////////////////////////////////////////
// version 01.02.2012
/*
// SIMPLEST USE
<script src="plot.js"></script> // my script
<canvas width="400" height="300" id="canv"></canvas> // created canvas element we draw on
<script> // configuration and initiation of plot
var myp = new MakeDraw(); // create drawing object
myp.id="canv"; // tell the object where would you like to draw
myp.data=data; // pass 1dimentional array with data you'd like to draw
myp.plot(); // drawing
</script>
*/
////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
options : 
this.horizontalNR = number -> number of horizontal lines ( myp.horizontalNR = number )
this.fSize=10; -> font size
this.textColor = 'rgba(red,green,blue,alpha)';
this.plotColor = 'rgba(red,green,blue,alpha)';
this.gridColor = 'rgba(red,green,blue,alpha)'; 
this.bgColor = 'rgba(red,green,blue,alpha)'; 
this.enumerateV = 1; -> vertical axis numeration
this.enumerateH = 1; -> horizontal axis numeration
this.lineWidthP = px -> plot line thickness
this.lineWidthG = px -> grid lines thickness
this.adjustTrimmer = nr -> adjusts auto precision 
*/
////////////////////////////////////////////////////////////////////////////////////////////////////////
function MakeDraw() {
///////////////////////// optional
this.prepSurface = function () {
document.write("<canvas id='canvas"+this.id+"' width="+this.width+" height="+this.height+"></canvas>");
}
/////////////////////////
//////// prepares user interface
this.prepUI = function () {
var canvas=document.getElementById(this.id);
var ctx=canvas.getContext('2d');
ctx.font = (this.fSize+"px sans-serif");
return ctx;
}
//////// determines spacing between horizontal and vertical lines
this.spacing = function (orientation,number) {
var canvas=document.getElementById(this.id);

if (orientation == "horizontal") {
var spac=((canvas.width-this.offsetL-this.offsetR)/number);
} else {
var spac=(canvas.height-(this.fSize*2))/number;
}

return spac;
}
//////// draws grid
this.drawGrid = function() {
var canvas=document.getElementById(this.id);
var hei = canvas.height-this.fSize*2;
var wid = canvas.width-this.offsetR;
var spacH,spacV;
var precalc;

ctx.fillStyle = this.bgColor;
ctx.fillRect(0,0,canvas.width,canvas.height);
ctx.lineWidth = this.lineWidthG;
ctx.strokeStyle = this.gridColor;
ctx.beginPath();


spacH = this.spacing("horizontal",this.data.length);
	for(var i=0, len = this.data.length; i<len;i++) {	
	precalc = i*spacH+this.offsetL;
	ctx.moveTo(precalc,1);
	ctx.lineTo(precalc,hei);
	}
	
spacV = this.spacing("vertical",this.horizontalNR);
	for(var i=0;i<this.horizontalNR;i++) {	
	precalc = i*spacV;
	ctx.moveTo(1,precalc);
	ctx.lineTo(wid-spacH,precalc);
	}

ctx.stroke(); 
ctx.closePath();
}
//////// stores minimal and maximal value in an array and returns it
this.getDataRange = function() {
var arr=new Array(0,0,0);
arr[0]=arr[1]=this.data[0];
	for(var i=1, len = this.data.length; i<len;i++) {
	if (this.data[i]<arr[0]) arr[0]=this.data[i];
	if (this.data[i]>arr[1]) arr[1]=this.data[i];
	}
return arr;
}
//////// modifies spacing basing on length of labels ( number )
function determineSpacing(number) {
var spacing=0;
	do {
	number/=10;
	spacing++;
	} while (number>1);
return spacing;
}
//////// determines offsets
this.determineOffsets = function() {
var range=this.getDataRange();
if (this.enumerateV) this.offsetL = (determineSpacing(this.dataTrimmer)+2+determineSpacing(Math.max(Math.abs(range[0]),Math.abs(range[1]))))*this.fSize*0.6;
if (this.enumerateP) this.offsetR = (determineSpacing(this.dataTrimmer)+2+determineSpacing(Math.abs(this.data[this.data.length-1])))*this.fSize*0.6;
if (this.enumerateH) this.offsetR = Math.max(this.offsetR,(determineSpacing(this.data.length))*this.fSize*0.6);
}
//////// draws linear graph and enumerates axes/curve
this.drawGraphLinear = function() { 
var canvas=document.getElementById(this.id);
var spacHoriz = this.spacing("horizontal",this.data.length);
var spacVertic = this.spacing("vertical",this.horizontalNR);
var hei=canvas.height-2*(spacVertic+this.fSize);

var range=this.getDataRange();
var totalRange=range[1]-range[0];
var verticalCoefficient=hei/totalRange;
var lookupTable = new Array();

	for(var i=0, len = this.data.length; i<len;i++) {
	lookupTable[i]=hei-(this.data[i]-range[0])*verticalCoefficient+spacVertic;
	}
	
ctx.lineWidth = this.lineWidthP;
ctx.strokeStyle = this.plotColor; 
ctx.beginPath();

	ctx.moveTo(this.offsetL,lookupTable[0]);
	for(var i=1, len = this.data.length; i<len;i++) {
	ctx.lineTo(i*spacHoriz+this.offsetL,lookupTable[i]);
	}
	
	

ctx.stroke(); 
ctx.closePath();
ctx.fillStyle = this.textColor;

	if (this.enumerateP) {
		for(var i=0, len = this.data.length; i<len;i++) {
			if (this.data[i]<this.data[i+1] && this.data[i-1]>this.data[i]) {
			ctx.fillText(Math.round(this.data[i]*this.dataTrimmer)/this.dataTrimmer,i*spacHoriz+this.offsetL,lookupTable[i]+12);
			} else {
			ctx.fillText(Math.round(this.data[i]*this.dataTrimmer)/this.dataTrimmer,i*spacHoriz+this.offsetL,lookupTable[i]-7);
			}
		}
	}
	if (this.enumerateH) {
	var spaceNeeded = this.data.length*this.fSize*determineSpacing(this.data.length);
		if (spaceNeeded < canvas.width) {
			for(var i=0, len = this.data.length; i<len;i++) {
			ctx.fillText(i+1,spacHoriz*i+this.offsetL-4,hei+2*spacVertic+10);
			}
		} else {
			for(var i=0; i<4;i++) {
			j=i*(this.data.length-1)/3;
			ctx.fillText(Math.round(j)+1,spacHoriz*j+this.offsetL-4,hei+2*spacVertic+10);
			}
		}
	}
	if (this.enumerateV) {
		for(var i=0; i<this.horizontalNR-1;i++) {
		var nrSpacing=totalRange/(this.horizontalNR-2);
		ctx.fillText(Math.round((range[1]-i*nrSpacing)*this.dataTrimmer)/this.dataTrimmer,1,(i+1)*spacVertic-this.fSize*0.5);
		}
	}
}
////////
this.trimData = function () {
	var range = this.getDataRange();
	var totalRange = range[1]-range[0];
	var spac = determineSpacing(100/totalRange); 
	this.dataTrimmer = Math.round(Math.pow(10,(spac+this.adjustTrimmer-1)));
}
//////// 
this.data;
this.horizontalNR = 10;
this.fSize=10;
this.lineWidthP=2;
this.lineWidthG=1;
this.dataTrimmer;
this.adjustTrimmer=0;
this.textColor = 'rgba(100,100,100,1)';
this.plotColor =  'rgba(200,100,100,1)';
this.gridColor =  'rgba(0,0,0,0.1)'; 
this.bgColor = 'rgba(255,255,255,1)'; 
this.enumerateV = 1;
this.enumerateH = 1;
this.enumerateP = 1;
this.offsetL=0;
this.offsetR=0;
this.id;
this.height;
this.width;
var ctx;
////////////////////// 
this.plot= function() {
//////
var canvas=document.getElementById(this.id);
this.width=canvas.width;
this.height=canvas.height;
//////
ctx=this.prepUI();
this.trimData();
this.determineOffsets();
this.drawGrid();
this.drawGraphLinear();
}
//////////////////////
}
////////////////////////////////////////