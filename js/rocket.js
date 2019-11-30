(function(){function r(e,n,t){function o(i,f){if(!n[i]){if(!e[i]){var c="function"==typeof require&&require;if(!f&&c)return c(i,!0);if(u)return u(i,!0);var a=new Error("Cannot find module '"+i+"'");throw a.code="MODULE_NOT_FOUND",a}var p=n[i]={exports:{}};e[i][0].call(p.exports,function(r){var n=e[i][1][r];return o(n||r)},p,p.exports,r,e,n,t)}return n[i].exports}for(var u="function"==typeof require&&require,i=0;i<t.length;i++)o(t[i]);return o}return r})()({1:[function(require,module,exports){

//import regression from 'regression'
var ode45 = require('ode45-cash-karp')
var linear = require('everpolate').linear
var eq, valid, hmax, peaked = false;

 // Original JavaScript code by Chirp Internet: www.chirp.com.au
  // Please acknowledge use of this code by including this header.

  window.requestAnimationFrame = window.requestAnimationFrame || window.mozRequestAnimationFrame || window.webkitRequestAnimationFrame || window.msRequestAnimationFrame;

  var field = document.getElementById("field");
  var ball = document.getElementById("ball");

  var maxX = field.clientWidth - ball.offsetWidth;
  var maxY = field.clientHeight - ball.offsetHeight;

  var duration = 20; // seconds
  var gridSize = 100; // pixels

  var start = null;

  var t = 0;
  var doneAnimating = false;


function getDesignSpan(){
    var designSpan = Number(document.getElementById('finSpan').value);
    return designSpan;
    // return 0.02;
 }
 function getDesignChord(){
    var designChord = Number(document.getElementById('finChord').value);
    return designChord;
    // return 0.2;
 }
 function getDesignLocation(){
    var designLocation = Number(document.getElementById('finLocation').value);
    return designLocation;
    // return 9.92;
 }
 function getWind(){
    var wind = Number(document.getElementById('wind').value);
    return wind;
    // return 200;
 }
 function getDesignFuselage(){
    var designFuselage = Number(document.getElementById('rocketLength').value);
    return designFuselage;
    // return 8.18;
 }
   function getRocketDiameter(){
    var rocketDiameter = Number(document.getElementById('rocketDiameter').value);
    return rocketDiameter;
    // return 0.038;
 }


window.simulate = function(){
    hmax = rocket_cost();
    if(valid)
      var outputText = document.getElementById('output_text').innerHTML = "Maximum height: " + Math.trunc(hmax) + " ft";
 }


// function rocket_ode(t,y){
    var rocketODE = function (dydt,y,t){


    var design_span = getDesignSpan();
    var design_chord = getDesignChord();
    var design_location = getDesignLocation();
    var wind = getWind();
    var design_fuselage = getDesignFuselage();
    var rocketDiameter = getRocketDiameter();

    // define constants

    var g = 9.8;                // gravity (assume constant for all altitudes)
    var rho = 1.225;            // density (assume constant for all altitudes)

    // change later
    var I = 2000;               // moment of inertia 

    var mm = 1.0;               // mass of motor
    var me = 0.5;               // mass of electronics bay
    var mc = 0.1;               // mass of parachute
    var mn = 0.25;              // mass of nosecone
    var mb = 0.0;               // mass of ballast
    var rhofus = 0.241;         // density (kg/m) of fuselage 
    var rhofin = 1800;          // density (kg/m^3) of fin

    // define coefficients of lift and drag 
    
    // // define coefficients of lift and drag
    // var CLofus = 0;       // CLo for fuselage
    // var CLofin = 0;       // CLo for fin

    // var CLafus =  0.125;       // CLalpha for fuselage
    // var CLafin = 0.125;     // CLalpha for fin

    // var CDofus = 0.015;     // CDo for fuselage
    // var CDofin = 0.015;    // CDo for fin

    // var CDafus = 0.0175;     // CDalpha for fuselage
    // var CDafin = 0.0175;  // CDalpha for fin

 var CLofus = 0;      
var CLofin = 0;      
var CLafus = 0;     
var CLafin = 0.1;   
var CDofus = 0.7;    
var CDofin = 0.015;   
var CDafin = 0.0001; 
var CDafus = 0.0;     
var CDafin = 0.0001; 


    // approximate thickness of airfoil
    var design_thickness = design_chord/6;


    // compute mass of fins and fuselage
    var mfus = rhofus*design_fuselage;
    var mfin = rhofin*(design_span*design_chord*design_thickness);

    // compute reference area for fuselage
    var A = Math.PI*Math.pow(rocketDiameter/2,2);

    // compute reference area for fin
    var S = design_span * design_chord;

    // compute body-axis velocities
    var u = y[2]+wind*Math.cos(y[4]);
    var w = y[3]+wind*Math.sin(y[4]);
    // console.log("W IS: "+w);
    // console.log("U IS: "+u);

    // compute total velocity
    var V = Math.sqrt(Math.pow(u,2) + Math.pow(w,2));

    // compute angle of attack 
    //   (set alpha=0 when forward velocity=0 to avoid divide-by-zero error);
    var alpha = 0;
    if (u != 0){
        alpha = Math.atan(w/u);
    }
    // console.log(" U IS " + u);

    // compute lift and drag at this velocity and angle-of-attack
    var L = 0.5*rho*Math.pow(V,2)*(A*(CLofus + CLafus*alpha) + S*(CLofin + CLafin*alpha));
    var D = 0.5*rho*Math.pow(V,2)*(A*(CDofus + CDafus*alpha) + S*(CDofin + CDafin*Math.pow(alpha,2)));

    // generate the thrust and mass of propellant
    var T  = rocket_thrust(t);
    // console.log(" T IS: " + t);
    var mp = rocket_mass(t);
    // console.log(" mp IS: " + mp);

    // generate weight
    var m = mp + mm + me + mc + mn + mfus + mfin;
    var W = m*g;

    // compute center of pressure of fins as distance from nose to quarter-chord
    var xcp = design_location + 0.25*design_chord;

    // compute center of gravity as distance from nose to cg
    var xcg = 4.5; // TODO: this will change


    // end simulation if we've hit ground already
    if (y[1] > 1){
        y.fill(0);
    }



    // generate state derivatives 
    dydt[0]= y[2]*Math.cos(y[4]) + y[3]*Math.sin(y[4]);
    dydt[1] = -1*y[2]*Math.sin(y[4]) + y[3]*Math.cos(y[4]);
    dydt[2] = (-1*W*Math.sin(y[4]) + T + L*Math.sin(alpha) - D*Math.cos(alpha))/m;
    dydt[3] = (W*Math.cos(y[4]) - L*Math.cos(alpha) - D*Math.sin(alpha))/m;
    dydt[4] = y[5];
    dydt[5] = -1*(xcp-xcg)*(L*Math.cos(alpha) + D*Math.sin(alpha))/I;
}

function rocket_cost(){

    // set the simulation length 
    var tfinal = 20;


    var dt0 = 1.5, y0 = [0, 0, 0, 0, Math.PI/2, 0];
    var integrator = ode45( y0, rocketODE, 0, dt0);

    // Integrate up to tmax:
    var t = [], y = [], count = 0, array = [], tempArr = [];
    while(integrator.step(tfinal)) {
      // Store the solution at this timestep:
      t.push( integrator.t );
      // Array.prototype.push.apply(tempArr, integrator.t);
      y.push( integrator.y );
      Array.prototype.push.apply(array, integrator.y);
    }


    var i, j, arr = [], tArr = [];

    while(array.length>0){
        arr.push(array.splice(0,6));
    }

    while(tempArr.length>0){
     tArr.push(tempArr.splice(0,6));
    }

// converting altitude/velocicty from meters to feet

    //console.log("LENGTH IS "+ y.length);
    for(i = 0; i<arr.length; i++){
        for(j = 0; j<4; j++){
            arr[i][j] = (arr[i][j]/0.3048);
        }
    }
    //console.log("y  4 4 " +y[4]);

    // Convert angles/rates from radians to degrees
    for(i = 0; i<arr.length; i++){
        for(j = 4; j<=5; j++){
            // console.log("DEBUGGING when i is " + i+ " "+y[i][j]);
            arr[i][j] = ((arr[i][j]*180)/Math.PI);
            // console.log("DEBUGGING when j is " + j+" "+y[i][j]);
        }
    }


    var h = [], dist = [];

    for(i = 0; i<arr.length; i++){
            h[i] = -1*arr[i][1];
            dist[i] = -1*arr[i][0];
    }

    // check maximum altitudes
    var hmax = Math.max(...h);

    var temp = [], data = [], foo = [], counter = 0;
    for(i = 0; i<h.length; i++){
        if(0<h[i]<1 && counter == 0){
           temp.push(t[i], h[i]);
           counter++;
        }
        if(h[i] == hmax){
          temp.push(t[i], h[i]);
          counter = 0;
        }
    }

    while(temp.length>0){
        data.push(temp.splice(0,2));
    }

    foo.push([0,0]);
    foo.push([7,882]);
    foo.push([16,0]);

    eq = regression.polynomial(data, { order: 2 });

    if(eq.string == "y = 0x^2 + 0x + 0"){
      document.getElementById("error-input").style.display = "block";
      valid=false;
    }
    else{
      document.getElementById("error-input").style.display = "none";
      valid=true;
    }
    // const equation = regression.polynomial([[0, 1], [32, 67], [12, 79]]);
    
     
    // draw(h, t);
    if(valid){
      step();
      return hmax;
    }

}


function rocket_thrust(t){

    // time data from GorillaMotor
    var tdata = [0, .015, .033, .049, .167, .310, .512, .618, .839, 0.899, 0.98, 1.02, 1.05, 1.26303, 1.3];

    // thrust data from GorillaMotors
    //
    var ThrustData = [0, 113.792, 193.101, 172.412, 256.893, 296.548, 303.445, 305.169, 298.272, 291.376, 266.893, 244.652, 222.411, 10.345, 0];

    var T = 0;
    if (t <Math.max(...tdata)){
        T = linear(t, tdata, ThrustData);
    }
    T = T-0;
    return T;
}




function rocket_mass(t){

// time data from GorillaMotor
    var tdata = [0, .015, .033, .049, .167, .310, .512, .618, .839, 0.899, 0.98, 1.02, 1.05, 1.26303, 1.3];

    // mass data from GorillaMotors
    //
mdata = [182/1000, 181.502/1000, 179.937/1000, 178.205/1000, 163.498/1000, 140.504/1000, 105.3/1000, 86.6684/1000, 47.9627/1000, 37.6479/1000, 24.5108/1000, 18.5712/1000, 14.5038/1000, 0.111027/1000, 0];

var mP = 0;
if (t <Math.max(...tdata)){
    //interpolation
    mP = linear(t, tdata, mdata)
    // console.log("INSIDE");
}

mP = mP - 0;

return mP;
}



function draw(t,yArr){

      // compile the expression once


      // evaluate the expression repeatedly for different values of x
      

      // render the plot using plotly
      const trace1 = {
        x: yArr,
        y: t,
        type: 'scatter'
      }
      const data = [trace1]
      Plotly.newPlot('plot', data)
 }



 

  var step  = function(timestamp){
    var progress, x, y, y2, t;

    if(start == null) {
      start = timestamp;
    }

    $('input').each(function() {
        if(!$(this).val()){
            $("#no-input").css("display", "block");
            $("retry-button").css("display", "none");
            $("simulate-button").css("display", "inline");
            valid = false;
           return false;
        }
        else{
            $("#no-input").css("display", "none");
            $("retry-button").css("display", "inline");
            $("simulate-button").css("display", "none");
            valid = true;
        }
    });

    if(valid){
    progress = (timestamp - start) / duration / 1000; // percent


    t = (timestamp - start)/1000;
    var time = document.getElementById('time').innerHTML = "Time: " + Math.trunc(t) + " s";

    x = (progress * maxX/gridSize)+1; // x = ƒ(t)

    y = ((eq.equation[0])*(Math.pow(t,2))+eq.equation[1]*t)*(1/500); // y = ƒ(x)
    yy = ((eq.equation[0])*(Math.pow(t,2))+eq.equation[1]*t);

    var parachute = document.getElementById("parachute");
    parachute.style.left =  Math.min(maxX, gridSize * x) -170+ "px";
    parachute.style.bottom = (gridSize * y)+4 + "px";

    ball.style.left =  Math.min(maxX, gridSize * x) + "px";
    ball.style.bottom = (gridSize * y) + "px";

    if(!peaked)
      degree = (yy/hmax)*10;
    else
      degree = (yy/hmax)+10*10;
    console.log(degree);

    ball.style.transform = "rotate("+degree+"deg)";
    if(yy>=hmax){
      peaked = true;
      parachute.style.visibility = "visible";
    }

    if(progress >= 1 || y<0) {
      doneAnimating = true;
      var x = document.getElementById("retry-button");
      var y = document.getElementById("simulate-button");
      x.style.display = "inline";
      y.style.display = "none";
    }

    if(!doneAnimating)
      requestAnimationFrame(step);
  }
}

window.retry  = function(){
  var x = document.getElementById("retry-button");
  var y = document.getElementById("simulate-button");
  // document.getElementById("input_form").reset();

  x.style.display = "none";
  y.style.display = "inline";
  ball.style.left =  Math.min(maxX, gridSize) + "px";
  ball.style.bottom = 0 + "px";
  ball.style.transform = "rotate("+0+"deg)";

  parachute.style.left =  Math.min(maxX, gridSize)-130 + "px";
  parachute.style.bottom = 0 + "px";
  parachute.style.visibility = "hidden";


  doneAnimating = false;
  start = null;
  peaked = false;

}




},{"everpolate":2,"ode45-cash-karp":8}],2:[function(require,module,exports){
'use strict';

// Expose module API:

module.exports.polynomial = require('./polynomial.js')
module.exports.linear = require('./linear.js')
module.exports.linearRegression = require('./linearRegression.js')
module.exports.step = require('./step.js')

},{"./linear.js":4,"./linearRegression.js":5,"./polynomial.js":6,"./step.js":7}],3:[function(require,module,exports){
'use strict';

/**
 * Makes argument to be an array if it's not
 *
 * @param input
 * @returns {Array}
 */

module.exports.makeItArrayIfItsNot = function (input) {
  return Object.prototype.toString.call( input ) !== '[object Array]'
    ? [input]
    : input
}

/**
 *
 * Utilizes bisection method to search an interval to which
 * point belongs to, then returns an index of left border
 * of the interval
 *
 * @param {Number} point
 * @param {Array} intervals
 * @returns {Number}
 */

module.exports.findIntervalLeftBorderIndex = function (point, intervals) {
  //If point is beyond given intervals
  if (point < intervals[0])
    return 0
  if (point > intervals[intervals.length - 1])
    return intervals.length - 1
  //If point is inside interval
  //Start searching on a full range of intervals
  var indexOfNumberToCompare 
    , leftBorderIndex = 0
    , rightBorderIndex = intervals.length - 1
  //Reduce searching range till it find an interval point belongs to using binary search
  while (rightBorderIndex - leftBorderIndex !== 1) {
    indexOfNumberToCompare = leftBorderIndex + Math.floor((rightBorderIndex - leftBorderIndex)/2)
    point >= intervals[indexOfNumberToCompare]
      ? leftBorderIndex = indexOfNumberToCompare
      : rightBorderIndex = indexOfNumberToCompare
  }
  return leftBorderIndex
}
},{}],4:[function(require,module,exports){
'use strict';

var help = require('./help')

module.exports = evaluateLinear

/**
 * Evaluates interpolating line/lines at the set of numbers
 * or at a single number for the function y=f(x)
 *
 * @param {Number|Array} pointsToEvaluate     number or set of numbers
 *                                            for which polynomial is calculated
 * @param {Array} functionValuesX             set of distinct x values
 * @param {Array} functionValuesY             set of distinct y=f(x) values
 * @returns {Array}
 */

function evaluateLinear (pointsToEvaluate, functionValuesX, functionValuesY) {
  var equations = []
  pointsToEvaluate = help.makeItArrayIfItsNot(pointsToEvaluate)
  pointsToEvaluate.forEach(function (point) {
    var index = help.findIntervalLeftBorderIndex(point, functionValuesX)
    if (index == functionValuesX.length - 1)
      index--
    equations.push(linearInterpolation(point, functionValuesX[index], functionValuesY[index]
      , functionValuesX[index + 1], functionValuesY[index + 1]))
  })
  return equations
}

/**
 *
 * Evaluates y-value at given x point for line that passes
 * through the points (x0,y0) and (y1,y1)
 *
 * @param x
 * @param x0
 * @param y0
 * @param x1
 * @param y1
 * @returns {Number}
 */

function linearInterpolation (x, x0, y0, x1, y1) {
  var a = (y1 - y0) / (x1 - x0)
  var b = -a * x0 + y0
  return a * x + b
}

},{"./help":3}],5:[function(require,module,exports){
'use strict';

module.exports = linearRegression

var help = require('./help')

/**
 * Computes Linear Regression slope, intercept, r-squared and returns
 * a function which can be used for evaluating linear regression
 * at a particular x-value
 *
 * @param functionValuesX {Array}
 * @param functionValuesY {Array}
 * @returns {Object}
 */

function linearRegression(functionValuesX, functionValuesY){
  var regression = {}
    , x = functionValuesX
    , y = functionValuesY
    , n = y.length
    , sum_x = 0
    , sum_y = 0
    , sum_xy = 0
    , sum_xx = 0
    , sum_yy = 0

  for (var i = 0; i < y.length; i++) {
    sum_x += x[i]
    sum_y += y[i]
    sum_xy += (x[i]*y[i])
    sum_xx += (x[i]*x[i])
    sum_yy += (y[i]*y[i])
  }

  regression.slope = (n * sum_xy - sum_x * sum_y) / (n*sum_xx - sum_x * sum_x)
  regression.intercept = (sum_y - regression.slope * sum_x)/n
  regression.rSquared = Math.pow((n*sum_xy - sum_x*sum_y)/Math.sqrt((n*sum_xx-sum_x*sum_x)*(n*sum_yy-sum_y*sum_y)),2)
  regression.evaluate = function (pointsToEvaluate) {
    var x = help.makeItArrayIfItsNot(pointsToEvaluate)
      , equation = []
      , that = this
    x.forEach(function (point) {
      equation.push(that.slope*point + that.intercept)
    })
    return equation
  }

  return regression
}

},{"./help":3}],6:[function(require,module,exports){
'use strict';

var help = require('./help')

module.exports = evaluatePolynomial

/**
 * Evaluates interpolating polynomial at the set of numbers
 * or at a single number for the function y=f(x)
 *
 * @param {Number|Array} pointsToEvaluate     number or set of numbers
 *                                            for which polynomial is calculated
 * @param {Array} functionValuesX             set of distinct x values
 * @param {Array} functionValuesY             set of distinct y=f(x) values
 * @returns {Array}                           interpolating polynomial
 */

 function evaluatePolynomial (pointsToEvaluate, functionValuesX, functionValuesY) {
  var equations = []
  pointsToEvaluate = help.makeItArrayIfItsNot(pointsToEvaluate)
  // evaluate the interpolating polynomial for each point
  pointsToEvaluate.forEach(function (point) {
    equations.push(nevillesIteratedInterpolation(point, functionValuesX, functionValuesY))
  })
  return equations
}

/**
 * Neville's Iterated Interpolation algorithm implementation
 * http://en.wikipedia.org/wiki/Neville's_algorithm <- for reference
 *
 * @param {Number} x                           number for which polynomial is calculated
 * @param {Array} X                            set of distinct x values
 * @param {Array} Y                            set of distinct y=f(x) values
 * @returns {number}                           interpolating polynomial
 */

function nevillesIteratedInterpolation (x, X, Y) {
  var Q = [Y]
  for (var i = 1; i < X.length; i++) {
    Q.push([])
    for (var j = 1; j <= i; j++) {
      Q[j][i] = ((x-X[i-j])*Q[j-1][i] - (x-X[i])*Q[j-1][i-1])/( X[i] - X[i-j] )
    }
  }
  return Q[j-1][i-1]
}

},{"./help":3}],7:[function(require,module,exports){
'use strict';

var help = require('./help')

module.exports = step

/**
 * Evaluates interpolating step function at the set of numbers
 * or at a single number
 *
 * @param {Number|Array} pointsToEvaluate     number or set of numbers
 *                                            for which step is calculated
 * @param {Array} functionValuesX             set of distinct x values
 * @param {Array} functionValuesY             set of distinct y=f(x) values
 * @returns {Array}
 */

function step (pointsToEvaluate, functionValuesX, functionValuesY) {
  return help.makeItArrayIfItsNot(pointsToEvaluate).map(function (point) {
    return functionValuesY[help.findIntervalLeftBorderIndex(point, functionValuesX)]
  })
}

},{"./help":3}],8:[function(require,module,exports){
'use strict'

module.exports = IntegratorFactory

function defaultErrorScaleFunction( i, dt, y, dydt ) {
  return Math.abs(y) + Math.abs(dt * dydt) + 1e-32
}

function defaultErrorReduceFunction( i, accumulatedError, errorEstimate ) {
  return Math.max( accumulatedError, Math.abs(errorEstimate))
}

function defaultErrorPostFunction( accumulatedError ) {
  return accumulatedError
}

function minMag (a, b) {
  return (a > 0 ? Math.min : Math.max)(a, b);
}

function maxMag (a, b) {
  return (a > 0 ? Math.max : Math.min)(a, b);
}

var Integrator = function Integrator( y0, deriv, t0, dt0, options ) {
  var opts = options || {}
  this.tol = opts.tol===undefined ? 1e-8 : opts.tol
  this.maxIncreaseFactor = opts.maxIncreaseFactor===undefined ? 10 : opts.maxIncreaseFactor
  this.maxDecreaseFactor = opts.maxDecreaseFactor===undefined ? 10 : opts.maxDecreaseFactor
  this.dtMinMag = opts.dtMinMag===undefined ? 0 : Math.abs(opts.dtMinMag)
  this.dtMaxMag = opts.dtMaxMag===undefined ? undefined : Math.abs(opts.dtMaxMag)
  this.verbose = opts.verbose===undefined ? true : !!opts.verbose;

  var logCnt = 0
  var maxLogs = 10
  var maxLogWarningIssued = false
  this.__log = function (method, msg) {
    if (!this.verbose) return;
    if (logCnt < maxLogs) {
      console.log('ode45-cash-karp::' + method + '(): ' + msg)
      logCnt++
    } else {
      if (!maxLogWarningIssued) {
        console.log('ode45-cash-karp: too many warnings. Silencing further output')
        maxLogWarningIssued = true
      }
    }
  }.bind(this)

  this.errorScaleFunction = opts.errorScaleFunction === undefined ? defaultErrorScaleFunction : opts.errorScaleFunction
  this.errorReduceFunction = opts.errorReduceFunction === undefined ? defaultErrorReduceFunction : opts.errorReduceFunction
  this.errorPostFunction = opts.errorPostFunction === undefined ? defaultErrorPostFunction : opts.errorPostFunction

  // This is technically a parameter, but I think the value of leaving this undocumented exceeds the
  // value of documenting this and only adding confusion. I can't imagine this will even need to be
  // modified.
  this.safetyFactor = opts.safetyFactor===undefined ? 0.9 : opts.safetyFactor

  // Bind variables to this:
  this.deriv = deriv
  this.y = y0
  this.n = this.y.length
  this.dt = dt0
  this.t = t0

  // Create a scratch array into which we compute the derivative:
  this._ctor = this.y.constructor

  this._errorScale = new this._ctor( this.n )
  this._w = new this._ctor( this.n )
  this._k1 = new this._ctor( this.n )
  this._k2 = new this._ctor( this.n )
  this._k3 = new this._ctor( this.n )
  this._k4 = new this._ctor( this.n )
  this._k5 = new this._ctor( this.n )
  this._k6 = new this._ctor( this.n )
}

Integrator.prototype._calculateK1 = function() {
  this.deriv( this._k1, this.y, this.t )

  return this
}

Integrator.prototype._calculateKs = function(dt) {
  var i

  //var a21 =  0.200000000000000000 // 1/5
  //var a31 =  0.075000000000000000 // 3/40
  //var a32 =  0.225000000000000000 // 9/40
  //var a41 =  0.300000000000000000 // 3/10
  //var a42 = -0.900000000000000000 // -9/10
  //var a43 =  1.200000000000000000 // 6/5
  //var a51 = -0.203703703703703703 // -11/54
  //var a52 =  2.500000000000000000 // 5/2
  //var a53 = -2.592592592592592592 // -70/27
  //var a54 =  1.296296296296296296 // 35/27
  //var a61 =  0.029495804398148148 // 1631/55296
  //var a62 =  0.341796875000000000 // 175/512
  //var a63 =  0.041594328703703703 // 575/13824
  //var a64 =  0.400345413773148148 // 44275/110592
  //var a65 =  0.061767578125000000 // 253/4096

  //var b1  =  0.000000000000000000 // 0
  //var b2  =  0.200000000000000000 // 1/5
  //var b3  =  0.300000000000000000 // 3/10
  //var b4  =  0.600000000000000000 // 3/5
  //var b5  =  1.000000000000000000 // 1
  //var b6  =  0.875000000000000000 // 7/8

  // Same for every step, so don't repeat:
  //this.deriv( this._k1, this.y, this.t )

  for(i=0; i<this.n; i++) {
    this._w[i] = this.y[i] + dt * (
      0.2 * this._k1[i] )
  }

  this.deriv( this._k2, this._w, this.t + dt * 0.2 )

  for(i=0; i<this.n; i++) {
    this._w[i] = this.y[i] + dt * (
      0.075 * this._k1[i] +
      0.225 * this._k2[i] )
  }

  this.deriv( this._k3, this._w, this.t + dt * 0.3 )

  for(i=0; i<this.n; i++) {
    this._w[i] = this.y[i] + dt * (
       0.3 * this._k1[i] +
      -0.9 * this._k2[i] +
       1.2 * this._k3[i] )
  }

  this.deriv( this._k4, this._w, this.t + dt * 0.6 )

  for(i=0; i<this.n; i++) {
    this._w[i] = this.y[i] + dt * (
      -0.203703703703703703 * this._k1[i] +
       2.5                  * this._k2[i] +
      -2.592592592592592592 * this._k3[i] +
       1.296296296296296296 * this._k4[i] )
  }

  this.deriv( this._k5, this._w, this.t + dt /* * b5 */ )

  for(i=0; i<this.n; i++) {
    this._w[i] = this.y[i] + dt * (
      0.029495804398148148 * this._k1[i] +
      0.341796875          * this._k2[i] +
      0.041594328703703703 * this._k3[i] +
      0.400345413773148148 * this._k4[i] +
      0.061767578125       * this._k5[i] )
  }

  this.deriv( this._k6, this._w, this.t + dt * 0.875 )

  return this
}

Integrator.prototype._calculateError = function(dt) {
  //var cs1 =  0.102177372685185185 // 2825/27648
  //var cs2 =  0.000000000000000000 // 0
  //var cs3 =  0.383907903439153439 // 18575/48384
  //var cs4 =  0.244592737268518518 // 13525/55296
  //var cs5 =  0.019321986607142857 // 277/14336
  //var cs6 =  0.250000000000000000 // 1/4

  //var dc1 =  0.004293774801587301 // cs1 - c1
  //var dc2 =  0.000000000000000000 // cs2 - c2
  //var dc3 = -0.018668586093857832 // cs3 - c3
  //var dc4 =  0.034155026830808080 // cs4 - c4
  //var dc5 =  0.019321986607142857 // cs5 - c5
  //var dc6 = -0.039102202145680406 // cs6 - c6

  var error = 0
  for(var i=0; i<this.n; i++) {
    error =  this.errorReduceFunction( i, error,
      dt * (
         0.004293774801587301 * this._k1[i] +
        -0.018668586093857832 * this._k3[i] +
         0.034155026830808080 * this._k4[i] +
         0.019321986607142857 * this._k5[i] +
        -0.039102202145680406 * this._k6[i]
      ) / this._errorScale[i]
    )
  }

  return this.errorPostFunction(error)
}

Integrator.prototype._update = function(dt) {
  //var c1  =  0.097883597883597883 // 37/378
  //var c2  =  0.000000000000000000 // 0
  //var c3  =  0.402576489533011272 // 250/621
  //var c4  =  0.210437710437710437 // 125/594
  //var c5  =  0.000000000000000000 // 0
  //var c6  =  0.289102202145680406 // 512/1771

  for(var i=0; i<this.n; i++) {
    this.y[i] += dt * (
      0.097883597883597883 * this._k1[i] +
      0.402576489533011272 * this._k3[i] +
      0.210437710437710437 * this._k4[i] +
      0.289102202145680406 * this._k6[i]
    )
  }
  this.t += dt
  return this
}

Integrator.prototype._calculateErrorScale = function(dt) {
  for(var i=0; i<this.n; i++) {
    this._errorScale[i] = this.errorScaleFunction(i, dt, this.y[i], this._k1[i])
  }
  return this
}

Integrator.prototype.step = function( tLimit ) {
  // Bail out early if we're *at* the limit:
  if (Math.abs(this.t - tLimit) < this.dt * 1e-10) {
    return false;
  }

  var thisDt = this.dt;

  // Don't integrate past a tLimit, if provided:
  if( tLimit !== undefined ) {
    thisDt = thisDt > 0 ? Math.min( tLimit - this.t, thisDt ) : Math.max( tLimit - this.t, thisDt )
  }

  // Limit the magnitude of dt to dtMaxMag
  if( this.dtMaxMag !== undefined && Math.abs( thisDt ) > this.dtMaxMag ) {
    this.__log('step', 'step greater than maximum stepsize requested. dt magnitude has been limited.')
    thisDt = thisDt > 0 ? this.dtMaxMag : -this.dtMaxMag
  }

  // Limit the magnitude of dt to dtMinMag
  if( this.dtMinMag !== undefined && Math.abs( thisDt ) < this.dtMinMag ) {
    this.__log('step', 'step smaller than minimum stepsize requested. dt magnitude has been limited.')
    thisDt = thisDt > 0 ? this.dtMinMag : -this.dtMinMag
  }

  // The first derivative doesn't change even if dt does, so only calculate this once:
  this._calculateK1()

  // The scale factor per-dimension probably doesn't need to change either across a single adaptive step:
  this._calculateErrorScale(thisDt)

  var error = Infinity
  var maxError = 0
  var nextDt
  var lowerDtLimitReached = false

  while(true) {

    // Calculate intermediate k's for the proposed step:
    this._calculateKs(thisDt)

    // Calculate the max error of the proposed step:
    error = this._calculateError(thisDt)

    if( error < this.tol || lowerDtLimitReached ) {
      // Success! Exit:
      break
    }

    if( ! Number.isFinite(error) ) {
      throw new Error('ode45-cash-karp::step() NaN encountered while integrating.')
    }

    // Failure. Adapt the timestep:
    nextDt = this.safetyFactor * thisDt * Math.pow( this.tol / error, 0.2 )

    // Cut the timestep, but not by more than maxDecreaseFactor
    thisDt = maxMag( thisDt / this.maxDecreaseFactor, nextDt )

    // If stepsize too small, finish off by taking the currently proposed step and logging a warning:
    if( this.dtMinMag !== undefined && Math.abs(thisDt) < this.dtMinMag ) {
      thisDt = this.dtMinMag * (thisDt > 0 ? 1 : -1);
      this.__log('step', 'minimum stepsize reached.')
      lowerDtLimitReached = true
    }
  }

  // Apply this update:
  this._update(thisDt)

  // Calculate the next timestep size:
  nextDt = this.safetyFactor * thisDt * Math.pow( this.tol / error, 0.25 )

  // Increase the timestep for the next time around, but not by more than the maxIncreaseFactor:
  this.dt = maxMag(this.dt / this.maxDecreaseFactor, minMag( this.dt * this.maxIncreaseFactor, nextDt ));

  if( tLimit !== undefined ) {
    return Math.abs(this.t - tLimit) > this.dt * 1e-8;
  } else {
    return true
  }
}

Integrator.prototype.steps = function( n, tLimit ) {
  for(var step=0; step<n; step++) {
    if( ! this.step(tLimit) ) return false;
  }
}

function IntegratorFactory( y0, deriv, t, dt, options ) {
  return new Integrator( y0, deriv, t, dt, options )
}

},{}]},{},[1]);
