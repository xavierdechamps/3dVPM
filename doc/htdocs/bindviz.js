// Note: For maximum-speed code, see "Optimizing Code" on the Emscripten wiki, https://github.com/kripken/emscripten/wiki/Optimizing-Code
// Note: Some Emscripten settings may limit the speed of the generated code.
try {
  this['Module'] = Module;
  Module.test;
} catch(e) {
  this['Module'] = Module = {};
}
// The environment setup code below is customized to use Module.
// *** Environment setup code ***
var ENVIRONMENT_IS_NODE = typeof process === 'object' && typeof require === 'function';
var ENVIRONMENT_IS_WEB = typeof window === 'object';
var ENVIRONMENT_IS_WORKER = typeof importScripts === 'function';
var ENVIRONMENT_IS_SHELL = !ENVIRONMENT_IS_WEB && !ENVIRONMENT_IS_NODE && !ENVIRONMENT_IS_WORKER;
if (typeof module === "object") {
  module.exports = Module;
}
if (ENVIRONMENT_IS_NODE) {
  // Expose functionality in the same simple way that the shells work
  // Note that we pollute the global namespace here, otherwise we break in node
  Module['print'] = function(x) {
    process['stdout'].write(x + '\n');
  };
  Module['printErr'] = function(x) {
    process['stderr'].write(x + '\n');
  };
  var nodeFS = require('fs');
  var nodePath = require('path');
  Module['read'] = function(filename, binary) {
    filename = nodePath['normalize'](filename);
    var ret = nodeFS['readFileSync'](filename);
    // The path is absolute if the normalized version is the same as the resolved.
    if (!ret && filename != nodePath['resolve'](filename)) {
      filename = path.join(__dirname, '..', 'src', filename);
      ret = nodeFS['readFileSync'](filename);
    }
    if (ret && !binary) ret = ret.toString();
    return ret;
  };
  Module['readBinary'] = function(filename) { return Module['read'](filename, true) };
  Module['load'] = function(f) {
    globalEval(read(f));
  };
  if (!Module['arguments']) {
    Module['arguments'] = process['argv'].slice(2);
  }
}
if (ENVIRONMENT_IS_SHELL) {
  Module['print'] = print;
  if (typeof printErr != 'undefined') Module['printErr'] = printErr; // not present in v8 or older sm
  Module['read'] = read;
  Module['readBinary'] = function(f) {
    return read(f, 'binary');
  };
  if (!Module['arguments']) {
    if (typeof scriptArgs != 'undefined') {
      Module['arguments'] = scriptArgs;
    } else if (typeof arguments != 'undefined') {
      Module['arguments'] = arguments;
    }
  }
}
if (ENVIRONMENT_IS_WEB && !ENVIRONMENT_IS_WORKER) {
  if (!Module['print']) {
    Module['print'] = function(x) {
      console.log(x);
    };
  }
  if (!Module['printErr']) {
    Module['printErr'] = function(x) {
      console.log(x);
    };
  }
}
if (ENVIRONMENT_IS_WEB || ENVIRONMENT_IS_WORKER) {
  Module['read'] = function(url) {
    var xhr = new XMLHttpRequest();
    xhr.open('GET', url, false);
    xhr.send(null);
    return xhr.responseText;
  };
  if (!Module['arguments']) {
    if (typeof arguments != 'undefined') {
      Module['arguments'] = arguments;
    }
  }
}
if (ENVIRONMENT_IS_WORKER) {
  // We can do very little here...
  var TRY_USE_DUMP = false;
  if (!Module['print']) {
    Module['print'] = (TRY_USE_DUMP && (typeof(dump) !== "undefined") ? (function(x) {
      dump(x);
    }) : (function(x) {
      // self.postMessage(x); // enable this if you want stdout to be sent as messages
    }));
  }
  Module['load'] = importScripts;
}
if (!ENVIRONMENT_IS_WORKER && !ENVIRONMENT_IS_WEB && !ENVIRONMENT_IS_NODE && !ENVIRONMENT_IS_SHELL) {
  // Unreachable because SHELL is dependant on the others
  throw 'Unknown runtime environment. Where are we?';
}
function globalEval(x) {
  eval.call(null, x);
}
if (!Module['load'] == 'undefined' && Module['read']) {
  Module['load'] = function(f) {
    globalEval(Module['read'](f));
  };
}
if (!Module['print']) {
  Module['print'] = function(){};
}
if (!Module['printErr']) {
  Module['printErr'] = Module['print'];
}
if (!Module['arguments']) {
  Module['arguments'] = [];
}
// *** Environment setup code ***
// Closure helpers
Module.print = Module['print'];
Module.printErr = Module['printErr'];
// Callbacks
if (!Module['preRun']) Module['preRun'] = [];
if (!Module['postRun']) Module['postRun'] = [];
// === Auto-generated preamble library stuff ===
//========================================
// Runtime code shared with compiler
//========================================
var Runtime = {
  stackSave: function () {
    return STACKTOP;
  },
  stackRestore: function (stackTop) {
    STACKTOP = stackTop;
  },
  forceAlign: function (target, quantum) {
    quantum = quantum || 4;
    if (quantum == 1) return target;
    if (isNumber(target) && isNumber(quantum)) {
      return Math.ceil(target/quantum)*quantum;
    } else if (isNumber(quantum) && isPowerOfTwo(quantum)) {
      var logg = log2(quantum);
      return '((((' +target + ')+' + (quantum-1) + ')>>' + logg + ')<<' + logg + ')';
    }
    return 'Math.ceil((' + target + ')/' + quantum + ')*' + quantum;
  },
  isNumberType: function (type) {
    return type in Runtime.INT_TYPES || type in Runtime.FLOAT_TYPES;
  },
  isPointerType: function isPointerType(type) {
  return type[type.length-1] == '*';
},
  isStructType: function isStructType(type) {
  if (isPointerType(type)) return false;
  if (isArrayType(type)) return true;
  if (/<?{ ?[^}]* ?}>?/.test(type)) return true; // { i32, i8 } etc. - anonymous struct types
  // See comment in isStructPointerType()
  return type[0] == '%';
},
  INT_TYPES: {"i1":0,"i8":0,"i16":0,"i32":0,"i64":0},
  FLOAT_TYPES: {"float":0,"double":0},
  or64: function (x, y) {
    var l = (x | 0) | (y | 0);
    var h = (Math.round(x / 4294967296) | Math.round(y / 4294967296)) * 4294967296;
    return l + h;
  },
  and64: function (x, y) {
    var l = (x | 0) & (y | 0);
    var h = (Math.round(x / 4294967296) & Math.round(y / 4294967296)) * 4294967296;
    return l + h;
  },
  xor64: function (x, y) {
    var l = (x | 0) ^ (y | 0);
    var h = (Math.round(x / 4294967296) ^ Math.round(y / 4294967296)) * 4294967296;
    return l + h;
  },
  getNativeTypeSize: function (type, quantumSize) {
    if (Runtime.QUANTUM_SIZE == 1) return 1;
    var size = {
      '%i1': 1,
      '%i8': 1,
      '%i16': 2,
      '%i32': 4,
      '%i64': 8,
      "%float": 4,
      "%double": 8
    }['%'+type]; // add '%' since float and double confuse Closure compiler as keys, and also spidermonkey as a compiler will remove 's from '_i8' etc
    if (!size) {
      if (type.charAt(type.length-1) == '*') {
        size = Runtime.QUANTUM_SIZE; // A pointer
      } else if (type[0] == 'i') {
        var bits = parseInt(type.substr(1));
        assert(bits % 8 == 0);
        size = bits/8;
      }
    }
    return size;
  },
  getNativeFieldSize: function (type) {
    return Math.max(Runtime.getNativeTypeSize(type), Runtime.QUANTUM_SIZE);
  },
  dedup: function dedup(items, ident) {
  var seen = {};
  if (ident) {
    return items.filter(function(item) {
      if (seen[item[ident]]) return false;
      seen[item[ident]] = true;
      return true;
    });
  } else {
    return items.filter(function(item) {
      if (seen[item]) return false;
      seen[item] = true;
      return true;
    });
  }
},
  set: function set() {
  var args = typeof arguments[0] === 'object' ? arguments[0] : arguments;
  var ret = {};
  for (var i = 0; i < args.length; i++) {
    ret[args[i]] = 0;
  }
  return ret;
},
  STACK_ALIGN: 8,
  getAlignSize: function (type, size, vararg) {
    // we align i64s and doubles on 64-bit boundaries, unlike x86
    if (type == 'i64' || type == 'double' || vararg) return 8;
    if (!type) return Math.min(size, 8); // align structures internally to 64 bits
    return Math.min(size || (type ? Runtime.getNativeFieldSize(type) : 0), Runtime.QUANTUM_SIZE);
  },
  calculateStructAlignment: function calculateStructAlignment(type) {
    type.flatSize = 0;
    type.alignSize = 0;
    var diffs = [];
    var prev = -1;
    type.flatIndexes = type.fields.map(function(field) {
      var size, alignSize;
      if (Runtime.isNumberType(field) || Runtime.isPointerType(field)) {
        size = Runtime.getNativeTypeSize(field); // pack char; char; in structs, also char[X]s.
        alignSize = Runtime.getAlignSize(field, size);
      } else if (Runtime.isStructType(field)) {
        size = Types.types[field].flatSize;
        alignSize = Runtime.getAlignSize(null, Types.types[field].alignSize);
      } else if (field[0] == 'b') {
        // bN, large number field, like a [N x i8]
        size = field.substr(1)|0;
        alignSize = 1;
      } else {
        throw 'Unclear type in struct: ' + field + ', in ' + type.name_ + ' :: ' + dump(Types.types[type.name_]);
      }
      if (type.packed) alignSize = 1;
      type.alignSize = Math.max(type.alignSize, alignSize);
      var curr = Runtime.alignMemory(type.flatSize, alignSize); // if necessary, place this on aligned memory
      type.flatSize = curr + size;
      if (prev >= 0) {
        diffs.push(curr-prev);
      }
      prev = curr;
      return curr;
    });
    type.flatSize = Runtime.alignMemory(type.flatSize, type.alignSize);
    if (diffs.length == 0) {
      type.flatFactor = type.flatSize;
    } else if (Runtime.dedup(diffs).length == 1) {
      type.flatFactor = diffs[0];
    }
    type.needsFlattening = (type.flatFactor != 1);
    return type.flatIndexes;
  },
  generateStructInfo: function (struct, typeName, offset) {
    var type, alignment;
    if (typeName) {
      offset = offset || 0;
      type = (typeof Types === 'undefined' ? Runtime.typeInfo : Types.types)[typeName];
      if (!type) return null;
      if (type.fields.length != struct.length) {
        printErr('Number of named fields must match the type for ' + typeName + ': possibly duplicate struct names. Cannot return structInfo');
        return null;
      }
      alignment = type.flatIndexes;
    } else {
      var type = { fields: struct.map(function(item) { return item[0] }) };
      alignment = Runtime.calculateStructAlignment(type);
    }
    var ret = {
      __size__: type.flatSize
    };
    if (typeName) {
      struct.forEach(function(item, i) {
        if (typeof item === 'string') {
          ret[item] = alignment[i] + offset;
        } else {
          // embedded struct
          var key;
          for (var k in item) key = k;
          ret[key] = Runtime.generateStructInfo(item[key], type.fields[i], alignment[i]);
        }
      });
    } else {
      struct.forEach(function(item, i) {
        ret[item[1]] = alignment[i];
      });
    }
    return ret;
  },
  dynCall: function (sig, ptr, args) {
    if (args && args.length) {
      return FUNCTION_TABLE[ptr].apply(null, args);
    } else {
      return FUNCTION_TABLE[ptr]();
    }
  },
  addFunction: function (func) {
    var table = FUNCTION_TABLE;
    var ret = table.length;
    table.push(func);
    table.push(0);
    return ret;
  },
  removeFunction: function (index) {
    var table = FUNCTION_TABLE;
    table[index] = null;
  },
  warnOnce: function (text) {
    if (!Runtime.warnOnce.shown) Runtime.warnOnce.shown = {};
    if (!Runtime.warnOnce.shown[text]) {
      Runtime.warnOnce.shown[text] = 1;
      Module.printErr(text);
    }
  },
  funcWrappers: {},
  getFuncWrapper: function (func, sig) {
    assert(sig);
    if (!Runtime.funcWrappers[func]) {
      Runtime.funcWrappers[func] = function() {
        return Runtime.dynCall(sig, func, arguments);
      };
    }
    return Runtime.funcWrappers[func];
  },
  UTF8Processor: function () {
    var buffer = [];
    var needed = 0;
    this.processCChar = function (code) {
      code = code & 0xff;
      if (needed) {
        buffer.push(code);
        needed--;
      }
      if (buffer.length == 0) {
        if (code < 128) return String.fromCharCode(code);
        buffer.push(code);
        if (code > 191 && code < 224) {
          needed = 1;
        } else {
          needed = 2;
        }
        return '';
      }
      if (needed > 0) return '';
      var c1 = buffer[0];
      var c2 = buffer[1];
      var c3 = buffer[2];
      var ret;
      if (c1 > 191 && c1 < 224) {
        ret = String.fromCharCode(((c1 & 31) << 6) | (c2 & 63));
      } else {
        ret = String.fromCharCode(((c1 & 15) << 12) | ((c2 & 63) << 6) | (c3 & 63));
      }
      buffer.length = 0;
      return ret;
    }
    this.processJSString = function(string) {
      string = unescape(encodeURIComponent(string));
      var ret = [];
      for (var i = 0; i < string.length; i++) {
        ret.push(string.charCodeAt(i));
      }
      return ret;
    }
  },
  stackAlloc: function (size) { var ret = STACKTOP;STACKTOP = (STACKTOP + size)|0;STACKTOP = ((((STACKTOP)+7)>>3)<<3); return ret; },
  staticAlloc: function (size) { var ret = STATICTOP;STATICTOP = (STATICTOP + size)|0;STATICTOP = ((((STATICTOP)+7)>>3)<<3); return ret; },
  dynamicAlloc: function (size) { var ret = DYNAMICTOP;DYNAMICTOP = (DYNAMICTOP + size)|0;DYNAMICTOP = ((((DYNAMICTOP)+7)>>3)<<3); if (DYNAMICTOP >= TOTAL_MEMORY) enlargeMemory();; return ret; },
  alignMemory: function (size,quantum) { var ret = size = Math.ceil((size)/(quantum ? quantum : 8))*(quantum ? quantum : 8); return ret; },
  makeBigInt: function (low,high,unsigned) { var ret = (unsigned ? (((low)>>>(0))+(((high)>>>(0))*4294967296)) : (((low)>>>(0))+(((high)|(0))*4294967296))); return ret; },
  GLOBAL_BASE: 8,
  QUANTUM_SIZE: 4,
  __dummy__: 0
}
//========================================
// Runtime essentials
//========================================
var __THREW__ = 0; // Used in checking for thrown exceptions.
var setjmpId = 1; // Used in setjmp/longjmp
var setjmpLabels = {};
var ABORT = false; // whether we are quitting the application. no code should run after this. set in exit() and abort()
var undef = 0;
// tempInt is used for 32-bit signed values or smaller. tempBigInt is used
// for 32-bit unsigned values or more than 32 bits. TODO: audit all uses of tempInt
var tempValue, tempInt, tempBigInt, tempInt2, tempBigInt2, tempPair, tempBigIntI, tempBigIntR, tempBigIntS, tempBigIntP, tempBigIntD;
var tempI64, tempI64b;
var tempRet0, tempRet1, tempRet2, tempRet3, tempRet4, tempRet5, tempRet6, tempRet7, tempRet8, tempRet9;
function abort(text) {
  Module.print(text + ':\n' + (new Error).stack);
  ABORT = true;
  throw "Assertion: " + text;
}
function assert(condition, text) {
  if (!condition) {
    abort('Assertion failed: ' + text);
  }
}
var globalScope = this;
// C calling interface. A convenient way to call C functions (in C files, or
// defined with extern "C").
//
// Note: LLVM optimizations can inline and remove functions, after which you will not be
//       able to call them. Closure can also do so. To avoid that, add your function to
//       the exports using something like
//
//         -s EXPORTED_FUNCTIONS='["_main", "_myfunc"]'
//
// @param ident      The name of the C function (note that C++ functions will be name-mangled - use extern "C")
// @param returnType The return type of the function, one of the JS types 'number', 'string' or 'array' (use 'number' for any C pointer, and
//                   'array' for JavaScript arrays and typed arrays).
// @param argTypes   An array of the types of arguments for the function (if there are no arguments, this can be ommitted). Types are as in returnType,
//                   except that 'array' is not possible (there is no way for us to know the length of the array)
// @param args       An array of the arguments to the function, as native JS values (as in returnType)
//                   Note that string arguments will be stored on the stack (the JS string will become a C string on the stack).
// @return           The return value, as a native JS value (as in returnType)
function ccall(ident, returnType, argTypes, args) {
  return ccallFunc(getCFunc(ident), returnType, argTypes, args);
}
Module["ccall"] = ccall;
// Returns the C function with a specified identifier (for C++, you need to do manual name mangling)
function getCFunc(ident) {
  try {
    var func = globalScope['Module']['_' + ident]; // closure exported function
    if (!func) func = eval('_' + ident); // explicit lookup
  } catch(e) {
  }
  assert(func, 'Cannot call unknown function ' + ident + ' (perhaps LLVM optimizations or closure removed it?)');
  return func;
}
// Internal function that does a C call using a function, not an identifier
function ccallFunc(func, returnType, argTypes, args) {
  var stack = 0;
  function toC(value, type) {
    if (type == 'string') {
      if (value === null || value === undefined || value === 0) return 0; // null string
      if (!stack) stack = Runtime.stackSave();
      var ret = Runtime.stackAlloc(value.length+1);
      writeStringToMemory(value, ret);
      return ret;
    } else if (type == 'array') {
      if (!stack) stack = Runtime.stackSave();
      var ret = Runtime.stackAlloc(value.length);
      writeArrayToMemory(value, ret);
      return ret;
    }
    return value;
  }
  function fromC(value, type) {
    if (type == 'string') {
      return Pointer_stringify(value);
    }
    assert(type != 'array');
    return value;
  }
  var i = 0;
  var cArgs = args ? args.map(function(arg) {
    return toC(arg, argTypes[i++]);
  }) : [];
  var ret = fromC(func.apply(null, cArgs), returnType);
  if (stack) Runtime.stackRestore(stack);
  return ret;
}
// Returns a native JS wrapper for a C function. This is similar to ccall, but
// returns a function you can call repeatedly in a normal way. For example:
//
//   var my_function = cwrap('my_c_function', 'number', ['number', 'number']);
//   alert(my_function(5, 22));
//   alert(my_function(99, 12));
//
function cwrap(ident, returnType, argTypes) {
  var func = getCFunc(ident);
  return function() {
    return ccallFunc(func, returnType, argTypes, Array.prototype.slice.call(arguments));
  }
}
Module["cwrap"] = cwrap;
// Sets a value in memory in a dynamic way at run-time. Uses the
// type data. This is the same as makeSetValue, except that
// makeSetValue is done at compile-time and generates the needed
// code then, whereas this function picks the right code at
// run-time.
// Note that setValue and getValue only do *aligned* writes and reads!
// Note that ccall uses JS types as for defining types, while setValue and
// getValue need LLVM types ('i8', 'i32') - this is a lower-level operation
function setValue(ptr, value, type, noSafe) {
  type = type || 'i8';
  if (type.charAt(type.length-1) === '*') type = 'i32'; // pointers are 32-bit
    switch(type) {
      case 'i1': HEAP8[(ptr)]=value; break;
      case 'i8': HEAP8[(ptr)]=value; break;
      case 'i16': HEAP16[((ptr)>>1)]=value; break;
      case 'i32': HEAP32[((ptr)>>2)]=value; break;
      case 'i64': (tempI64 = [value>>>0,Math.min(Math.floor((value)/4294967296), 4294967295)>>>0],HEAP32[((ptr)>>2)]=tempI64[0],HEAP32[(((ptr)+(4))>>2)]=tempI64[1]); break;
      case 'float': HEAPF32[((ptr)>>2)]=value; break;
      case 'double': HEAPF64[((ptr)>>3)]=value; break;
      default: abort('invalid type for setValue: ' + type);
    }
}
Module['setValue'] = setValue;
// Parallel to setValue.
function getValue(ptr, type, noSafe) {
  type = type || 'i8';
  if (type.charAt(type.length-1) === '*') type = 'i32'; // pointers are 32-bit
    switch(type) {
      case 'i1': return HEAP8[(ptr)];
      case 'i8': return HEAP8[(ptr)];
      case 'i16': return HEAP16[((ptr)>>1)];
      case 'i32': return HEAP32[((ptr)>>2)];
      case 'i64': return HEAP32[((ptr)>>2)];
      case 'float': return HEAPF32[((ptr)>>2)];
      case 'double': return HEAPF64[((ptr)>>3)];
      default: abort('invalid type for setValue: ' + type);
    }
  return null;
}
Module['getValue'] = getValue;
var ALLOC_NORMAL = 0; // Tries to use _malloc()
var ALLOC_STACK = 1; // Lives for the duration of the current function call
var ALLOC_STATIC = 2; // Cannot be freed
var ALLOC_DYNAMIC = 3; // Cannot be freed except through sbrk
var ALLOC_NONE = 4; // Do not allocate
Module['ALLOC_NORMAL'] = ALLOC_NORMAL;
Module['ALLOC_STACK'] = ALLOC_STACK;
Module['ALLOC_STATIC'] = ALLOC_STATIC;
Module['ALLOC_DYNAMIC'] = ALLOC_DYNAMIC;
Module['ALLOC_NONE'] = ALLOC_NONE;
// allocate(): This is for internal use. You can use it yourself as well, but the interface
//             is a little tricky (see docs right below). The reason is that it is optimized
//             for multiple syntaxes to save space in generated code. So you should
//             normally not use allocate(), and instead allocate memory using _malloc(),
//             initialize it with setValue(), and so forth.
// @slab: An array of data, or a number. If a number, then the size of the block to allocate,
//        in *bytes* (note that this is sometimes confusing: the next parameter does not
//        affect this!)
// @types: Either an array of types, one for each byte (or 0 if no type at that position),
//         or a single type which is used for the entire block. This only matters if there
//         is initial data - if @slab is a number, then this does not matter at all and is
//         ignored.
// @allocator: How to allocate memory, see ALLOC_*
function allocate(slab, types, allocator, ptr) {
  var zeroinit, size;
  if (typeof slab === 'number') {
    zeroinit = true;
    size = slab;
  } else {
    zeroinit = false;
    size = slab.length;
  }
  var singleType = typeof types === 'string' ? types : null;
  var ret;
  if (allocator == ALLOC_NONE) {
    ret = ptr;
  } else {
    ret = [_malloc, Runtime.stackAlloc, Runtime.staticAlloc, Runtime.dynamicAlloc][allocator === undefined ? ALLOC_STATIC : allocator](Math.max(size, singleType ? 1 : types.length));
  }
  if (zeroinit) {
    var ptr = ret, stop;
    assert((ret & 3) == 0);
    stop = ret + (size & ~3);
    for (; ptr < stop; ptr += 4) {
      HEAP32[((ptr)>>2)]=0;
    }
    stop = ret + size;
    while (ptr < stop) {
      HEAP8[((ptr++)|0)]=0;
    }
    return ret;
  }
  if (singleType === 'i8') {
    if (slab.subarray || slab.slice) {
      HEAPU8.set(slab, ret);
    } else {
      HEAPU8.set(new Uint8Array(slab), ret);
    }
    return ret;
  }
  var i = 0, type, typeSize, previousType;
  while (i < size) {
    var curr = slab[i];
    if (typeof curr === 'function') {
      curr = Runtime.getFunctionIndex(curr);
    }
    type = singleType || types[i];
    if (type === 0) {
      i++;
      continue;
    }
    if (type == 'i64') type = 'i32'; // special case: we have one i32 here, and one i32 later
    setValue(ret+i, curr, type);
    // no need to look up size unless type changes, so cache it
    if (previousType !== type) {
      typeSize = Runtime.getNativeTypeSize(type);
      previousType = type;
    }
    i += typeSize;
  }
  return ret;
}
Module['allocate'] = allocate;
function Pointer_stringify(ptr, /* optional */ length) {
  // Find the length, and check for UTF while doing so
  var hasUtf = false;
  var t;
  var i = 0;
  while (1) {
    t = HEAPU8[(((ptr)+(i))|0)];
    if (t >= 128) hasUtf = true;
    else if (t == 0 && !length) break;
    i++;
    if (length && i == length) break;
  }
  if (!length) length = i;
  var ret = '';
  if (!hasUtf) {
    var MAX_CHUNK = 1024; // split up into chunks, because .apply on a huge string can overflow the stack
    var curr;
    while (length > 0) {
      curr = String.fromCharCode.apply(String, HEAPU8.subarray(ptr, ptr + Math.min(length, MAX_CHUNK)));
      ret = ret ? ret + curr : curr;
      ptr += MAX_CHUNK;
      length -= MAX_CHUNK;
    }
    return ret;
  }
  var utf8 = new Runtime.UTF8Processor();
  for (i = 0; i < length; i++) {
    t = HEAPU8[(((ptr)+(i))|0)];
    ret += utf8.processCChar(t);
  }
  return ret;
}
Module['Pointer_stringify'] = Pointer_stringify;
// Memory management
var PAGE_SIZE = 4096;
function alignMemoryPage(x) {
  return ((x+4095)>>12)<<12;
}
var HEAP;
var HEAP8, HEAPU8, HEAP16, HEAPU16, HEAP32, HEAPU32, HEAPF32, HEAPF64;
var STATIC_BASE = 0, STATICTOP = 0, staticSealed = false; // static area
var STACK_BASE = 0, STACKTOP = 0, STACK_MAX = 0; // stack area
var DYNAMIC_BASE = 0, DYNAMICTOP = 0; // dynamic area handled by sbrk
function enlargeMemory() {
  // TOTAL_MEMORY is the current size of the actual array, and DYNAMICTOP is the new top.
  while (TOTAL_MEMORY <= DYNAMICTOP) { // Simple heuristic. Override enlargeMemory() if your program has something more optimal for it
    TOTAL_MEMORY = alignMemoryPage(2*TOTAL_MEMORY);
  }
  assert(TOTAL_MEMORY <= Math.pow(2, 30)); // 2^30==1GB is a practical maximum - 2^31 is already close to possible negative numbers etc.
  var oldHEAP8 = HEAP8;
  var buffer = new ArrayBuffer(TOTAL_MEMORY);
  Module['HEAP8'] = HEAP8 = new Int8Array(buffer);
  Module['HEAP16'] = HEAP16 = new Int16Array(buffer);
  Module['HEAP32'] = HEAP32 = new Int32Array(buffer);
  Module['HEAPU8'] = HEAPU8 = new Uint8Array(buffer);
  Module['HEAPU16'] = HEAPU16 = new Uint16Array(buffer);
  Module['HEAPU32'] = HEAPU32 = new Uint32Array(buffer);
  Module['HEAPF32'] = HEAPF32 = new Float32Array(buffer);
  Module['HEAPF64'] = HEAPF64 = new Float64Array(buffer);
  HEAP8.set(oldHEAP8);
}
var TOTAL_STACK = Module['TOTAL_STACK'] || 5242880;
var TOTAL_MEMORY = Module['TOTAL_MEMORY'] || 16777216;
var FAST_MEMORY = Module['FAST_MEMORY'] || 2097152;
// Initialize the runtime's memory
// check for full engine support (use string 'subarray' to avoid closure compiler confusion)
assert(!!Int32Array && !!Float64Array && !!(new Int32Array(1)['subarray']) && !!(new Int32Array(1)['set']),
       'Cannot fallback to non-typed array case: Code is too specialized');
var buffer = new ArrayBuffer(TOTAL_MEMORY);
HEAP8 = new Int8Array(buffer);
HEAP16 = new Int16Array(buffer);
HEAP32 = new Int32Array(buffer);
HEAPU8 = new Uint8Array(buffer);
HEAPU16 = new Uint16Array(buffer);
HEAPU32 = new Uint32Array(buffer);
HEAPF32 = new Float32Array(buffer);
HEAPF64 = new Float64Array(buffer);
// Endianness check (note: assumes compiler arch was little-endian)
HEAP32[0] = 255;
assert(HEAPU8[0] === 255 && HEAPU8[3] === 0, 'Typed arrays 2 must be run on a little-endian system');
Module['HEAP'] = HEAP;
Module['HEAP8'] = HEAP8;
Module['HEAP16'] = HEAP16;
Module['HEAP32'] = HEAP32;
Module['HEAPU8'] = HEAPU8;
Module['HEAPU16'] = HEAPU16;
Module['HEAPU32'] = HEAPU32;
Module['HEAPF32'] = HEAPF32;
Module['HEAPF64'] = HEAPF64;
function callRuntimeCallbacks(callbacks) {
  while(callbacks.length > 0) {
    var callback = callbacks.shift();
    if (typeof callback == 'function') {
      callback();
      continue;
    }
    var func = callback.func;
    if (typeof func === 'number') {
      if (callback.arg === undefined) {
        Runtime.dynCall('v', func);
      } else {
        Runtime.dynCall('vi', func, [callback.arg]);
      }
    } else {
      func(callback.arg === undefined ? null : callback.arg);
    }
  }
}
var __ATINIT__ = []; // functions called during startup
var __ATMAIN__ = []; // functions called when main() is to be run
var __ATEXIT__ = []; // functions called during shutdown
var runtimeInitialized = false;
function ensureInitRuntime() {
  if (runtimeInitialized) return;
  runtimeInitialized = true;
  callRuntimeCallbacks(__ATINIT__);
}
function preMain() {
  callRuntimeCallbacks(__ATMAIN__);
}
function exitRuntime() {
  callRuntimeCallbacks(__ATEXIT__);
}
// Tools
// This processes a JS string into a C-line array of numbers, 0-terminated.
// For LLVM-originating strings, see parser.js:parseLLVMString function
function intArrayFromString(stringy, dontAddNull, length /* optional */) {
  var ret = (new Runtime.UTF8Processor()).processJSString(stringy);
  if (length) {
    ret.length = length;
  }
  if (!dontAddNull) {
    ret.push(0);
  }
  return ret;
}
Module['intArrayFromString'] = intArrayFromString;
function intArrayToString(array) {
  var ret = [];
  for (var i = 0; i < array.length; i++) {
    var chr = array[i];
    if (chr > 0xFF) {
      chr &= 0xFF;
    }
    ret.push(String.fromCharCode(chr));
  }
  return ret.join('');
}
Module['intArrayToString'] = intArrayToString;
// Write a Javascript array to somewhere in the heap
function writeStringToMemory(string, buffer, dontAddNull) {
  var array = intArrayFromString(string, dontAddNull);
  var i = 0;
  while (i < array.length) {
    var chr = array[i];
    HEAP8[(((buffer)+(i))|0)]=chr
    i = i + 1;
  }
}
Module['writeStringToMemory'] = writeStringToMemory;
function writeArrayToMemory(array, buffer) {
  for (var i = 0; i < array.length; i++) {
    HEAP8[(((buffer)+(i))|0)]=array[i];
  }
}
Module['writeArrayToMemory'] = writeArrayToMemory;
function unSign(value, bits, ignore, sig) {
  if (value >= 0) {
    return value;
  }
  return bits <= 32 ? 2*Math.abs(1 << (bits-1)) + value // Need some trickery, since if bits == 32, we are right at the limit of the bits JS uses in bitshifts
                    : Math.pow(2, bits)         + value;
}
function reSign(value, bits, ignore, sig) {
  if (value <= 0) {
    return value;
  }
  var half = bits <= 32 ? Math.abs(1 << (bits-1)) // abs is needed if bits == 32
                        : Math.pow(2, bits-1);
  if (value >= half && (bits <= 32 || value > half)) { // for huge values, we can hit the precision limit and always get true here. so don't do that
                                                       // but, in general there is no perfect solution here. With 64-bit ints, we get rounding and errors
                                                       // TODO: In i64 mode 1, resign the two parts separately and safely
    value = -2*half + value; // Cannot bitshift half, as it may be at the limit of the bits JS uses in bitshifts
  }
  return value;
}
if (!Math['imul']) Math['imul'] = function(a, b) {
  var ah  = a >>> 16;
  var al = a & 0xffff;
  var bh  = b >>> 16;
  var bl = b & 0xffff;
  return (al*bl + ((ah*bl + al*bh) << 16))|0;
};
// A counter of dependencies for calling run(). If we need to
// do asynchronous work before running, increment this and
// decrement it. Incrementing must happen in a place like
// PRE_RUN_ADDITIONS (used by emcc to add file preloading).
// Note that you can add dependencies in preRun, even though
// it happens right before run - run will be postponed until
// the dependencies are met.
var runDependencies = 0;
var runDependencyTracking = {};
var calledInit = false, calledRun = false;
var runDependencyWatcher = null;
function addRunDependency(id) {
  runDependencies++;
  if (Module['monitorRunDependencies']) {
    Module['monitorRunDependencies'](runDependencies);
  }
  if (id) {
    assert(!runDependencyTracking[id]);
    runDependencyTracking[id] = 1;
  } else {
    Module.printErr('warning: run dependency added without ID');
  }
}
Module['addRunDependency'] = addRunDependency;
function removeRunDependency(id) {
  runDependencies--;
  if (Module['monitorRunDependencies']) {
    Module['monitorRunDependencies'](runDependencies);
  }
  if (id) {
    assert(runDependencyTracking[id]);
    delete runDependencyTracking[id];
  } else {
    Module.printErr('warning: run dependency removed without ID');
  }
  if (runDependencies == 0) {
    if (runDependencyWatcher !== null) {
      clearInterval(runDependencyWatcher);
      runDependencyWatcher = null;
    } 
    // If run has never been called, and we should call run (INVOKE_RUN is true, and Module.noInitialRun is not false)
    if (!calledRun && shouldRunNow) run();
  }
}
Module['removeRunDependency'] = removeRunDependency;
Module["preloadedImages"] = {}; // maps url to image data
Module["preloadedAudios"] = {}; // maps url to audio data
function addPreRun(func) {
  if (!Module['preRun']) Module['preRun'] = [];
  else if (typeof Module['preRun'] == 'function') Module['preRun'] = [Module['preRun']];
  Module['preRun'].push(func);
}
var awaitingMemoryInitializer = false;
function loadMemoryInitializer(filename) {
  function applyData(data) {
    HEAPU8.set(data, STATIC_BASE);
    runPostSets();
  }
  // always do this asynchronously, to keep shell and web as similar as possible
  addPreRun(function() {
    if (ENVIRONMENT_IS_NODE || ENVIRONMENT_IS_SHELL) {
      applyData(Module['readBinary'](filename));
    } else {
      Browser.asyncLoad(filename, function(data) {
        applyData(data);
      }, function(data) {
        throw 'could not load memory initializer ' + filename;
      });
    }
  });
  awaitingMemoryInitializer = false;
}
// === Body ===
STATIC_BASE = 8;
STATICTOP = STATIC_BASE + 194032;
var _stdout;
var _stdin;
var _stderr;
__ATINIT__ = __ATINIT__.concat([
  { func: function() { __GLOBAL__I_a() } },
  { func: function() { __GLOBAL__I_a5924() } },
  { func: function() { __GLOBAL__I_a6109() } }
]);
var ___fsmu8;
var ___dso_handle;
var __ZTVN10__cxxabiv120__si_class_type_infoE;
var __ZTVN10__cxxabiv119__pointer_type_infoE;
var __ZTVN10__cxxabiv117__class_type_infoE;
var __ZTIt;
var __ZTIs;
var __ZTIm;
var __ZTIl;
var __ZTIj;
var __ZTIi;
var __ZTIh;
var __ZTIf;
var __ZTId;
var __ZTIc;
var __ZTIa;
var _stdout = _stdout=allocate([0,0,0,0,0,0,0,0], "i8", ALLOC_STATIC);
var _stdin = _stdin=allocate([0,0,0,0,0,0,0,0], "i8", ALLOC_STATIC);
var _stderr = _stderr=allocate([0,0,0,0,0,0,0,0], "i8", ALLOC_STATIC);
__ZTVN10__cxxabiv120__si_class_type_infoE=allocate([0,0,0,0,136,225,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], "i8", ALLOC_STATIC);
__ZTVN10__cxxabiv119__pointer_type_infoE=allocate([0,0,0,0,152,225,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], "i8", ALLOC_STATIC);
__ZTVN10__cxxabiv117__class_type_infoE=allocate([0,0,0,0,184,225,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], "i8", ALLOC_STATIC);
__ZTIt=allocate([0,0,0,0,0,0,0,0], "i8", ALLOC_STATIC);
__ZTIs=allocate([0,0,0,0,0,0,0,0], "i8", ALLOC_STATIC);
__ZTIm=allocate([0,0,0,0,0,0,0,0], "i8", ALLOC_STATIC);
__ZTIl=allocate([0,0,0,0,0,0,0,0], "i8", ALLOC_STATIC);
__ZTIj=allocate([0,0,0,0,0,0,0,0], "i8", ALLOC_STATIC);
__ZTIi=allocate([0,0,0,0,0,0,0,0], "i8", ALLOC_STATIC);
__ZTIh=allocate([0,0,0,0,0,0,0,0], "i8", ALLOC_STATIC);
__ZTIf=allocate([0,0,0,0,0,0,0,0], "i8", ALLOC_STATIC);
__ZTId=allocate([0,0,0,0,0,0,0,0], "i8", ALLOC_STATIC);
__ZTIc=allocate([0,0,0,0,0,0,0,0], "i8", ALLOC_STATIC);
__ZTIa=allocate([0,0,0,0,0,0,0,0], "i8", ALLOC_STATIC);
/* memory initializer */ allocate([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,32,1,0,0,238,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,58,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,28,1,0,0,0,0,0,0,134,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,70,1,0,0,0,0,0,0,90,2,0,0,18,0,0,0,86,1,0,0,76,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,240,0,0,0,0,1,0,0,16,1,0,0,32,1,0,0,48,1,0,0,64,1,0,0,80,1,0,0,96,1,0,0,0,1,0,0,0,1,0,0,64,1,0,0,64,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,248,63,51,51,51,51,51,51,211,63,0,0,0,0,0,0,0,0,60,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,64,164,2,0,152,149,2,0,200,135,2,0,240,120,2,0,112,106,2,0,104,94,2,0,8,81,2,0,208,65,2,0,136,54,2,0,160,46,2,0,216,38,2,0,128,32,2,0,208,24,2,0,0,17,2,0,192,6,2,0,104,2,2,0,80,1,0,0,0,0,0,0,44,3,0,0,200,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,36,4,0,0,118,3,0,0,0,0,0,0,0,0,0,0,50,1,0,0,0,0,0,0,124,0,0,0,242,4,0,0,122,4,0,0,230,1,0,0,70,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,204,2,0,0,126,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,216,2,0,0,0,0,0,0,138,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,238,1,0,0,0,0,0,0,194,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,176,3,0,0,0,0,0,0,144,2,0,0,244,3,0,0,196,2,0,0,144,0,0,0,138,4,0,0,0,0,0,0,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,88,168,53,205,59,78,213,63,37,117,2,154,8,27,218,63,0,0,0,0,0,0,224,63,0,0,0,0,0,0,224,63,44,212,154,230,29,167,234,63,106,222,113,138,142,228,232,63,88,168,53,205,59,78,213,63,88,168,53,205,59,78,213,63,88,168,53,205,59,78,213,63,0,0,0,0,0,0,224,63,93,220,70,3,120,11,226,63,0,0,0,0,0,0,208,63,88,168,53,205,59,78,213,63,0,0,0,0,0,0,208,63,211,188,227,20,29,201,209,63,0,0,0,0,0,0,224,63,0,0,0,0,0,0,224,63,0,0,0,0,0,0,224,63,0,0,0,0,0,0,224,63,0,0,0,0,0,0,224,63,0,0,0,0,0,0,224,63,0,0,0,0,0,0,224,63,0,0,0,0,0,0,224,63,0,0,0,0,0,0,224,63,0,0,0,0,0,0,224,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,93,220,70,3,120,11,226,63,93,220,70,3,120,11,226,63,93,220,70,3,120,11,226,63,13,113,172,139,219,104,220,63,100,93,220,70,3,120,237,63,210,111,95,7,206,25,231,63,16,122,54,171,62,87,229,63,16,122,54,171,62,87,229,63,210,111,95,7,206,25,231,63,120,11,36,40,126,140,227,63,181,21,251,203,238,201,225,63,210,111,95,7,206,25,231,63,210,111,95,7,206,25,231,63,88,168,53,205,59,78,213,63,136,133,90,211,188,227,216,63,210,111,95,7,206,25,231,63,120,11,36,40,126,140,227,63,196,66,173,105,222,113,236,63,210,111,95,7,206,25,231,63,210,111,95,7,206,25,231,63,181,21,251,203,238,201,225,63,210,111,95,7,206,25,231,63,16,122,54,171,62,87,229,63,181,21,251,203,238,201,225,63,120,11,36,40,126,140,227,63,210,111,95,7,206,25,231,63,210,111,95,7,206,25,231,63,134,56,214,197,109,52,238,63,210,111,95,7,206,25,231,63,210,111,95,7,206,25,231,63,120,11,36,40,126,140,227,63,88,168,53,205,59,78,213,63,211,188,227,20,29,201,209,63,88,168,53,205,59,78,213,63,166,10,70,37,117,2,222,63,0,0,0,0,0,0,224,63,88,168,53,205,59,78,213,63,13,113,172,139,219,104,220,63,0,0,0,0,0,0,224,63,13,113,172,139,219,104,220,63,0,0,0,0,0,0,224,63,13,113,172,139,219,104,220,63,88,168,53,205,59,78,213,63,0,0,0,0,0,0,224,63,0,0,0,0,0,0,224,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,0,0,0,0,0,0,224,63,211,188,227,20,29,201,209,63,106,222,113,138,142,228,232,63,0,0,0,0,0,0,224,63,0,0,0,0,0,0,224,63,0,0,0,0,0,0,224,63,0,0,0,0,0,0,224,63,88,168,53,205,59,78,213,63,136,133,90,211,188,227,216,63,211,188,227,20,29,201,209,63,0,0,0,0,0,0,224,63,0,0,0,0,0,0,224,63,210,111,95,7,206,25,231,63,0,0,0,0,0,0,224,63,0,0,0,0,0,0,224,63,13,113,172,139,219,104,220,63,244,108,86,125,174,182,222,63,17,54,60,189,82,150,201,63,244,108,86,125,174,182,222,63,59,1,77,132,13,79,225,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,88,168,53,205,59,78,213,63,0,0,0,0,0,0,224,63,0,0,0,0,0,0,224,63,62,232,217,172,250,92,197,63,0,0,0,0,0,0,224,63,0,0,0,0,0,0,224,63,0,0,0,0,0,0,224,63,0,0,0,0,0,0,224,63,130,115,70,148,246,6,199,63,13,113,172,139,219,104,220,63,0,0,0,0,0,0,224,63,88,168,53,205,59,78,213,63,88,168,53,205,59,78,213,63,181,21,251,203,238,201,225,63,181,21,251,203,238,201,225,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,224,63,0,0,0,0,0,0,224,63,0,0,0,0,0,0,224,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,7,240,22,72,80,252,220,63,162,180,55,248,194,100,214,63,88,168,53,205,59,78,213,63,13,113,172,139,219,104,220,63,13,113,172,139,219,104,220,63,0,0,0,0,0,0,224,63,0,0,0,0,0,0,240,63,0,0,0,0,0,0,240,63,0,0,0,0,0,0,208,63,13,113,172,139,219,104,220,63,0,0,0,0,0,0,208,63,88,168,53,205,59,78,213,63,88,168,53,205,59,78,213,63,88,168,53,205,59,78,213,63,88,168,53,205,59,78,213,63,88,168,53,205,59,78,213,63,88,168,53,205,59,78,213,63,88,168,53,205,59,78,213,63,88,168,53,205,59,78,213,63,0,0,0,0,0,0,208,63,88,168,53,205,59,78,213,63,88,168,53,205,59,78,213,63,0,0,0,0,0,0,208,63,88,168,53,205,59,78,213,63,88,168,53,205,59,78,213,63,88,168,53,205,59,78,213,63,0,0,0,0,0,0,240,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,196,66,173,105,222,113,236,63,0,0,0,0,0,0,208,63,127,217,61,121,88,168,209,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,120,11,36,40,126,140,227,63,210,111,95,7,206,25,231,63,196,66,173,105,222,113,236,63,19,242,65,207,102,213,211,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,16,122,54,171,62,87,229,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,211,188,227,20,29,201,209,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,211,188,227,20,29,201,209,63,0,0,0,0,0,0,224,63,210,111,95,7,206,25,231,63,0,0,0,0,0,0,224,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,0,0,0,0,0,0,208,63,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,2,0,0,0,2,0,0,0,1,0,0,0,2,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,32,3,2,0,56,234,1,0,64,164,2,0,216,189,1,0,160,169,1,0,72,150,1,0,120,132,1,0,152,149,2,0,224,164,2,0,200,135,2,0,56,136,2,0,168,121,2,0,8,107,2,0,208,94,2,0,168,81,2,0,144,69,2,0,192,54,2,0,8,47,2,0,80,39,2,0,216,32,2,0,24,25,2,0,88,17,2,0,8,7,2,0,136,2,2,0,0,0,2,0,16,254,1,0,40,252,1,0,248,249,1,0,80,247,1,0,16,245,1,0,32,243,1,0,48,240,1,0,16,236,1,0,208,233,1,0,168,231,1,0,112,229,1,0,184,226,1,0,112,224,1,0,72,222,1,0,72,220,1,0,152,218,1,0,152,216,1,0,96,212,1,0,208,209,1,0,152,207,1,0,128,205,1,0,216,203,1,0,104,202,1,0,240,120,2,0,232,198,1,0,48,197,1,0,240,193,1,0,136,191,1,0,112,106,2,0,104,94,2,0,248,185,1,0,64,184,1,0,144,182,1,0,184,180,1,0,184,178,1,0,240,176,1,0,0,175,1,0,216,171,1,0,56,169,1,0,32,167,1,0,160,165,1,0,136,163,1,0,176,161,1,0,232,159,1,0,56,158,1,0,112,156,1,0,80,154,1,0,192,151,1,0,56,149,1,0,104,147,1,0,136,145,1,0,40,144,1,0,192,142,1,0,104,141,1,0,240,139,1,0,88,138,1,0,8,136,1,0,8,81,2,0,32,132,1,0,136,130,1,0,16,129,1,0,208,65,2,0,32,126,1,0,160,124,1,0,200,122,1,0,104,121,1,0,160,119,1,0,168,117,1,0,240,115,1,0,56,114,1,0,0,113,1,0,248,111,1,0,240,110,1,0,208,172,2,0,168,171,2,0,96,170,2,0,136,54,2,0,88,166,2,0,160,46,2,0,16,163,2,0,168,161,2,0,112,160,2,0,32,159,2,0,176,157,2,0,128,156,2,0,48,155,2,0,176,153,2,0,240,151,2,0,8,150,2,0,120,148,2,0,64,147,2,0,8,146,2,0,176,144,2,0,216,38,2,0,128,32,2,0,184,140,2,0,136,139,2,0,176,137,2,0,8,136,2,0,208,134,2,0,112,133,2,0,184,131,2,0,160,130,2,0,208,24,2,0,248,127,2,0,104,126,2,0,160,124,2,0,208,122,2,0,80,121,2,0,8,120,2,0,176,118,2,0,160,116,2,0,0,17,2,0,240,113,2,0,200,112,2,0,104,111,2,0,144,109,2,0,96,108,2,0,192,6,2,0,96,105,2,0,104,2,2,0,240,102,2,0,0,0,0,0,0,0,0,0,0,0,0,0,78,0,0,0,0,0,0,0,210,4,0,0,36,5,0,0,48,1,0,0,190,1,0,0,6,2,0,0,154,0,0,0,10,1,0,0,30,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,226,1,0,0,88,3,0,0,234,1,0,0,42,4,0,0,80,4,0,0,152,5,0,0,0,0,0,0,0,0,0,0,2,3,0,0,0,0,0,0,180,2,0,0,150,5,0,0,156,3,0,0,96,2,0,0,164,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,255,255,255,255,0,0,0,0,252,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,0,0,0,255,255,255,255,0,0,0,0,0,0,0,0,204,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,125,0,0,0,0,0,0,0,10,0,0,0,0,0,0,0,105,108,108,45,99,111,110,100,105,116,105,111,110,101,100,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,252,2,0,0,6,0,0,0,0,0,0,0,0,0,0,0,98,5,0,0,156,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,48,0,0,0,38,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,32,193,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,176,193,0,0,0,0,0,0,0,0,0,0,0,0,0,160,1,0,0,16,0,0,0,1,0,0,0,0,0,0,0,0,16,0,2,0,0,0,0,0,0,0,0,0,0,16,64,0,0,0,0,0,0,0,0,5,0,0,0,0,0,0,0,0,176,193,0,0,0,0,0,0,0,0,0,0,0,16,64,136,11,0,0,147,0,0,0,1,0,0,0,0,0,0,0,0,32,3,2,0,0,0,0,0,0,0,0,0,0,16,64,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,64,32,62,3,0,0,0,0,0,0,0,0,0,0,16,64,96,21,0,0,122,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,16,64,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,144,195,0,0,0,0,0,0,0,0,0,0,0,16,64,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,16,0,0,0,0,0,0,0,0,0,0,0,0,16,64,216,67,0,0,8,0,0,0,1,0,0,0,0,0,0,0,0,32,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5,0,0,0,0,0,0,0,168,1,0,0,244,0,0,0,68,1,0,0,150,2,0,0,6,4,0,0,28,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,24,49,2,0,224,96,2,0,200,83,2,0,128,72,2,0,88,56,2,0,0,0,0,0,1,0,0,0,2,0,0,0,3,0,0,0,4,0,0,0,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,130,3,0,0,24,3,0,0,42,0,0,0,0,0,0,0,172,5,0,0,0,0,0,0,82,3,0,0,142,0,0,0,24,5,0,0,188,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,170,0,0,0,206,2,0,0,242,2,0,0,160,0,0,0,212,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,168,0,0,0,0,0,0,0,134,0,0,0,102,0,0,0,60,0,0,0,12,4,0,0,232,2,0,0,28,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,248,176,1,0,8,175,1,0,224,171,1,0,96,65,2,0,48,167,1,0,176,165,1,0,152,163,1,0,192,161,1,0,96,65,2,0,248,159,1,0,72,158,1,0,96,65,2,0,136,156,1,0,96,154,1,0,208,151,1,0,72,149,1,0,120,147,1,0,152,145,1,0,56,144,1,0,208,142,1,0,120,141,1,0,0,140,1,0,104,138,1,0,24,136,1,0,248,133,1,0,48,132,1,0,144,130,1,0,24,129,1,0,168,127,1,0,56,126,1,0,176,124,1,0,216,122,1,0,120,121,1,0,176,119,1,0,96,65,2,0,184,117,1,0,8,116,1,0,72,114,1,0,16,113,1,0,96,65,2,0,8,112,1,0,0,111,1,0,224,172,2,0,184,171,2,0,176,119,1,0,96,65,2,0,112,170,2,0,160,168,2,0,96,166,2,0,128,164,2,0,32,163,2,0,176,161,2,0,128,160,2,0,40,159,2,0,192,157,2,0,144,156,2,0,64,155,2,0,96,65,2,0,192,153,2,0,0,152,2,0,24,150,2,0,128,148,2,0,72,147,2,0,96,65,2,0,16,146,2,0,192,144,2,0,32,143,2,0,240,141,2,0,200,140,2,0,152,139,2,0,192,137,2,0,16,136,2,0,224,134,2,0,128,133,2,0,200,131,2,0,168,130,2,0,176,119,1,0,96,65,2,0,112,129,2,0,0,128,2,0,120,126,2,0,208,142,1,0,96,65,2,0,176,124,2,0,224,122,2,0,88,121,2,0,24,120,2,0,192,118,2,0,168,116,2,0,56,115,2,0,248,113,2,0,208,112,2,0,120,111,2,0,208,142,1,0,96,65,2,0,152,109,2,0,104,108,2,0,208,106,2,0,112,105,2,0,24,104,2,0,0,103,2,0,248,101,2,0,232,100,2,0,176,119,1,0,96,65,2,0,16,100,2,0,168,98,2,0,80,97,2,0,32,96,2,0,160,94,2,0,72,93,2,0,8,92,2,0,232,90,2,0,208,89,2,0,112,88,2,0,88,87,2,0,176,119,1,0,96,65,2,0,232,85,2,0,40,84,2,0,96,65,2,0,248,82,2,0,80,81,2,0,24,80,2,0,40,79,2,0,80,78,2,0,72,77,2,0,96,76,2,0,120,75,2,0,120,74,2,0,96,65,2,0,80,73,2,0,96,65,2,0,248,70,2,0,16,66,2,0,192,64,2,0,216,63,2,0,200,62,2,0,160,61,2,0,176,119,1,0,96,65,2,0,64,60,2,0,96,65,2,0,64,59,2,0,24,58,2,0,184,56,2,0,136,55,2,0,168,54,2,0,48,54,2,0,192,53,2,0,208,142,1,0,96,65,2,0,240,52,2,0,96,65,2,0,64,52,2,0,176,51,2,0,32,51,2,0,80,50,2,0,88,49,2,0,40,48,2,0,216,46,2,0,96,65,2,0,208,45,2,0,40,45,2,0,104,44,2,0,168,43,2,0,16,43,2,0,120,42,2,0,224,41,2,0,32,41,2,0,96,65,2,0,40,40,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,36,64,0,0,0,0,0,0,89,64,0,0,0,0,0,136,195,64,0,0,0,0,132,215,151,65,0,128,224,55,121,195,65,67,23,110,5,181,181,184,147,70,245,249,63,233,3,79,56,77,50,29,48,249,72,119,130,90,60,191,115,127,221,79,21,117,216,189,1,0,184,190,1,0,152,149,2,0,200,135,2,0,56,136,2,0,216,116,1,0,88,165,2,0,8,151,2,0,160,136,2,0,168,121,2,0,208,94,2,0,80,95,2,0,48,82,2,0,24,70,2,0,192,54,2,0,8,47,2,0,24,25,2,0,24,33,2,0,16,254,1,0,16,245,1,0,48,240,1,0,176,2,2,0,112,229,1,0,184,226,1,0,112,224,1,0,64,250,1,0,72,222,1,0,56,245,1,0,64,243,1,0,72,240,1,0,120,236,1,0,128,205,1,0,192,231,1,0,104,202,1,0,240,193,1,0,136,191,1,0,104,222,1,0,112,220,1,0,192,218,1,0,208,216,1,0,136,212,1,0,232,209,1,0,184,207,1,0,200,205,1,0,0,204,1,0,160,202,1,0,192,200,1,0,88,199,1,0,112,197,1,0,8,194,1,0,168,191,1,0,120,189,1,0,104,187,1,0,32,186,1,0,104,184,1,0,104,94,2,0,208,180,1,0,248,185,1,0,8,177,1,0,184,178,1,0,216,171,1,0,176,161,1,0,64,167,1,0,88,138,1,0,184,163,1,0,32,132,1,0,16,129,1,0,104,158,1,0,208,65,2,0,32,126,1,0,160,124,1,0,128,149,1,0,144,147,1,0,200,122,1,0,160,119,1,0,168,117,1,0,240,115,1,0,56,114,1,0,0,113,1,0,64,136,1,0,24,134,1,0,88,132,1,0,248,111,1,0,136,54,2,0,192,127,1,0,88,126,1,0,208,124,1,0,0,123,1,0,136,121,1,0,192,119,1,0,168,161,2,0,112,160,2,0,32,159,2,0,128,156,2,0,64,147,2,0,8,146,2,0,16,173,2,0,128,32,2,0,128,170,2,0,8,136,2,0,160,166,2,0,112,133,2,0,120,163,2,0,160,130,2,0,208,24,2,0,248,127,2,0,104,126,2,0,168,156,2,0,8,120,2,0,176,118,2,0,40,152,2,0,160,116,2,0,240,113,2,0,104,111,2,0,48,146,2,0,216,144,2,0,144,109,2,0,0,142,2,0,96,108,2,0,192,6,2,0,104,2,2,0,240,102,2,0,188,1,0,0,0,0,0,0,196,1,0,0,138,5,0,0,72,2,0,0,208,3,0,0,186,4,0,0,4,1,0,0,96,0,0,0,120,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,54,0,0,0,166,3,0,0,118,2,0,0,140,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,62,2,0,0,0,0,0,0,114,5,0,0,22,2,0,0,112,4,0,0,50,4,0,0,210,3,0,0,0,0,0,0,112,47,2,0,0,121,1,0,16,170,2,0,0,0,0,0,0,0,0,0,4,0,0,0,200,154,2,0,0,0,0,0,0,0,0,0,16,55,2,0,0,121,1,0,16,170,2,0,0,0,0,0,0,126,2,0,5,0,0,0,200,154,2,0,0,0,0,0,80,111,2,0,64,70,2,0,0,121,1,0,184,85,2,0,0,0,0,0,0,0,0,0,6,0,0,0,200,154,2,0,120,70,2,0,0,0,0,0,168,39,2,0,0,121,1,0,184,85,2,0,0,0,0,0,0,126,2,0,7,0,0,0,200,154,2,0,120,70,2,0,80,111,2,0,240,233,1,0,176,41,2,0,184,85,2,0,0,0,0,0,0,0,0,0,10,0,0,0,88,35,2,0,120,70,2,0,0,0,0,0,224,226,1,0,176,41,2,0,184,85,2,0,0,0,0,0,80,111,2,0,11,0,0,0,88,35,2,0,120,70,2,0,80,111,2,0,136,229,1,0,176,41,2,0,216,10,2,0,0,0,0,0,0,0,0,0,8,0,0,0,88,35,2,0,0,0,0,0,0,0,0,0,200,231,1,0,176,41,2,0,216,10,2,0,0,0,0,0,80,111,2,0,9,0,0,0,88,35,2,0,0,0,0,0,80,111,2,0,144,7,2,0,144,7,2,0,0,0,0,0,0,0,0,0,0,0,0,0,12,0,0,0,232,254,1,0,0,0,0,0,0,0,0,0,200,17,2,0,144,7,2,0,120,70,2,0,0,0,0,0,0,0,0,0,14,0,0,0,232,254,1,0,120,70,2,0,0,0,0,0,208,2,2,0,144,7,2,0,120,70,2,0,0,0,0,0,0,126,2,0,15,0,0,0,232,254,1,0,120,70,2,0,80,111,2,0,104,248,1,0,144,7,2,0,0,0,0,0,0,0,0,0,0,126,2,0,13,0,0,0,232,254,1,0,0,0,0,0,80,111,2,0,64,0,2,0,64,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,16,0,0,0,200,154,2,0,0,0,0,0,0,0,0,0,64,254,1,0,64,0,2,0,120,70,2,0,0,0,0,0,0,0,0,0,18,0,0,0,200,154,2,0,120,70,2,0,0,0,0,0,88,250,1,0,64,0,2,0,120,70,2,0,0,0,0,0,0,126,2,0,19,0,0,0,200,154,2,0,120,70,2,0,80,111,2,0,128,243,1,0,64,0,2,0,0,0,0,0,200,234,1,0,0,0,0,0,20,0,0,0,200,154,2,0,0,0,0,0,0,0,0,0,136,247,1,0,64,0,2,0,120,70,2,0,200,234,1,0,0,0,0,0,22,0,0,0,200,154,2,0,120,70,2,0,0,0,0,0,104,240,1,0,64,0,2,0,120,70,2,0,200,234,1,0,0,126,2,0,23,0,0,0,200,154,2,0,120,70,2,0,80,111,2,0,88,245,1,0,64,0,2,0,0,0,0,0,200,234,1,0,0,126,2,0,21,0,0,0,200,154,2,0,0,0,0,0,80,111,2,0,104,252,1,0,64,0,2,0,0,0,0,0,0,0,0,0,0,126,2,0,17,0,0,0,200,154,2,0,0,0,0,0,80,111,2,0,160,224,1,0,88,221,1,0,120,70,2,0,0,0,0,0,0,0,0,0,26,0,0,0,88,35,2,0,120,70,2,0,0,0,0,0,232,218,1,0,88,221,1,0,120,70,2,0,0,0,0,0,80,111,2,0,27,0,0,0,88,35,2,0,120,70,2,0,80,111,2,0,112,222,1,0,88,221,1,0,0,0,0,0,0,0,0,0,80,111,2,0,25,0,0,0,88,35,2,0,0,0,0,0,80,111,2,0,120,220,1,0,88,221,1,0,232,210,1,0,0,0,0,0,0,0,0,0,24,0,0,0,88,35,2,0,0,0,0,0,0,0,0,0,144,212,1,0,152,206,1,0,120,70,2,0,0,0,0,0,0,0,0,0,30,0,0,0,88,35,2,0,120,70,2,0,0,0,0,0,192,207,1,0,152,206,1,0,120,70,2,0,0,0,0,0,80,111,2,0,31,0,0,0,88,35,2,0,120,70,2,0,80,111,2,0,240,209,1,0,152,206,1,0,0,0,0,0,0,0,0,0,80,111,2,0,29,0,0,0,88,35,2,0,0,0,0,0,80,111,2,0,216,216,1,0,152,206,1,0,232,210,1,0,0,0,0,0,0,0,0,0,28,0,0,0,88,35,2,0,0,0,0,0,0,0,0,0,8,204,1,0,8,204,1,0,0,0,0,0,0,0,0,0,0,0,0,0,32,0,0,0,32,198,1,0,0,0,0,0,0,0,0,0,64,33,2,0,96,192,1,0,120,70,2,0,0,0,0,0,0,0,0,0,2,0,0,0,88,35,2,0,120,70,2,0,0,0,0,0,128,25,2,0,96,192,1,0,120,70,2,0,0,0,0,0,80,111,2,0,3,0,0,0,88,35,2,0,120,70,2,0,80,111,2,0,152,236,1,0,96,192,1,0,0,0,0,0,0,0,0,0,80,111,2,0,1,0,0,0,88,35,2,0,0,0,0,0,80,111,2,0,208,205,1,0,96,192,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,88,35,2,0,0,0,0,0,0,0,0,0,24,185,1,0,120,183,1,0,184,181,1,0,0,0,0,0,80,111,2,0,33,0,0,0,88,35,2,0,0,0,0,0,80,111,2,0,168,202,1,0,240,177,1,0,0,0,0,0,0,0,0,0,0,0,0,0,34,0,0,0,32,198,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,98,2,0,0,232,3,0,0,162,3,0,0,116,5,0,0,204,1,0,0,108,5,0,0,184,126,1,0,24,125,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,50,0,0,0,232,3,0,0,162,3,0,0,106,1,0,0,0,0,0,0,82,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,100,111,116,32,112,105,99,32,112,108,117,103,105,110,58,32,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,166,1,0,0,164,0,0,0,0,0,0,0,0,0,0,0,62,1,0,0,124,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,102,1,0,0,0,0,0,0,158,4,0,0,228,2,0,0,54,3,0,0,12,3,0,0,92,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,240,63,0,0,0,0,0,0,240,63,0,0,0,0,0,0,240,63,0,0,0,0,0,0,240,63,0,0,0,0,0,0,240,63,0,0,0,0,0,0,240,63,0,0,0,0,0,0,240,63,0,0,0,0,0,0,240,63,0,0,0,0,0,0,240,63,0,0,0,0,0,0,240,63,0,0,0,0,0,0,240,63,0,0,0,0,0,0,240,63,184,5,0,0,64,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,3,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,154,153,153,153,153,153,217,191,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,13,0,0,0,0,0,0,0,0,1,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,0,0,0,0,0,0,0,0,1,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,1,0,0,0,1,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,16,0,0,0,0,0,0,0,0,1,0,0,0,7,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,22,0,0,0,0,0,0,0,0,1,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,19,0,0,0,0,0,0,0,0,1,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,18,0,0,0,0,0,0,0,0,1,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,11,0,0,0,0,0,0,0,0,1,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,23,0,0,0,0,0,0,0,0,1,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,21,0,0,0,0,0,0,0,0,1,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,20,0,0,0,0,0,0,0,0,1,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6,0,0,0,0,0,0,0,0,1,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0].concat([0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,51,51,51,51,51,51,227,63,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,14,0,0,0,0,0,0,0,0,1,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,25,0,0,0,0,0,0,0,0,1,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,24,0,0,0,0,0,0,0,0,1,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,128,102,64,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,128,102,64,154,153,153,153,153,153,217,191,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,5,0,0,0,0,0,0,0,0,0,0,0,0,128,102,64,123,20,174,71,225,122,228,191,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,17,0,0,0,0,0,0,0,0,1,0,0,0,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,123,20,174,71,225,122,228,191,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,0,0,0,0,0,0,0,0,1,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,12,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,51,51,51,51,51,51,211,191,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,2,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,128,70,64,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5,0,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7,0,0,0,0,0,0,0,0,1,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,1,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,15,0,0,0,0,1,0,0,0,1,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,128,70,64,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,24,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,24,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,82,1,0,0,240,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,162,2,0,0,88,4,0,0,46,5,0,0,88,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,0,0,0,4,0,0,0,0,0,0,0,88,2,0,0,170,5,0,0,88,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,58,2,0,0,88,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,132,0,0,0,0,0,0,0,0,0,0,0,192,1,0,0,0,0,0,0,0,0,0,0,194,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,45,57,57,57,57,57,57,57,57,57,57,57,57,57,57,57,46,57,57,0,0,0,0,0,0,0,0,0,0,0,0,0,154,153,153,153,153,153,169,63,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,216,0,0,0,2,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,94,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,0,0,0,4,0,0,0,0,0,0,0,234,4,0,0,70,5,0,0,24,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,162,5,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,24,252,1,0,8,0,0,0,3,0,0,0,240,249,1,0,64,247,1,0,11,0,0,0,6,0,0,0,160,197,1,0,24,243,1,0,2,0,0,0,1,0,0,0,40,240,1,0,8,236,1,0,4,0,0,0,2,0,0,0,200,233,1,0,160,231,1,0,4,0,0,0,4,0,0,0,104,229,1,0,176,226,1,0,5,0,0,0,5,0,0,0,104,224,1,0,64,222,1,0,4,0,0,0,7,0,0,0,64,220,1,0,144,218,1,0,5,0,0,0,9,0,0,0,144,216,1,0,88,212,1,0,4,0,0,0,10,0,0,0,200,209,1,0,1,208,209,210,211,212,213,214,215,216,217,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,0,0,0,8,0,0,0,0,0,0,0,170,3,0,0,170,2,0,0,108,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,16,1,0,0,48,1,0,0,176,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,208,131,1,0,1,0,0,0,224,1,0,0,192,15,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,216,186,1,0,1,0,0,0,208,2,0,0,224,15,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,80,191,1,0,1,0,0,0,224,13,0,0,0,16,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,160,197,1,0,1,0,0,0,136,17,0,0,32,16,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,48,4,2,0,1,0,0,0,72,23,0,0,64,16,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,64,202,1,0,255,255,255,255,200,33,0,0,96,16,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,40,203,1,0,1,0,0,0,104,48,0,0,128,16,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,144,203,1,0,1,0,0,0,248,67,0,0,160,16,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,80,235,1,0,1,0,0,0,96,83,0,0,192,16,0,0,4,0,0,0,112,211,1,0,1,0,0,0,64,0,0,0,160,15,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,184,205,1,0,112,93,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,96,206,1,0,128,93,1,0,160,204,1,0,208,92,1,0,15,0,0,0,120,205,1,0,1,0,0,0,152,82,0,0,0,0,0,0,16,0,0,0,16,104,2,0,1,0,0,0,152,82,0,0,0,0,0,0,17,0,0,0,88,5,2,0,1,0,0,0,152,82,0,0,0,0,0,0,17,0,0,0,104,235,1,0,1,0,0,0,152,82,0,0,0,0,0,0,17,0,0,0,136,211,1,0,1,0,0,0,152,82,0,0,0,0,0,0,19,0,0,0,16,191,1,0,1,0,0,0,184,82,0,0,0,0,0,0,20,0,0,0,200,170,1,0,1,0,0,0,184,82,0,0,0,0,0,0,21,0,0,0,56,151,1,0,1,0,0,0,184,82,0,0,0,0,0,0,21,0,0,0,136,133,1,0,1,0,0,0,184,82,0,0,0,0,0,0,21,0,0,0,64,117,1,0,1,0,0,0,184,82,0,0,0,0,0,0,22,0,0,0,200,165,2,0,1,0,0,0,136,82,0,0,0,0,0,0,23,0,0,0,152,151,2,0,1,0,0,0,136,82,0,0,0,0,0,0,24,0,0,0,96,137,2,0,1,0,0,0,136,82,0,0,0,0,0,0,24,0,0,0,120,122,2,0,1,0,0,0,136,82,0,0,0,0,0,0,24,0,0,0,32,108,2,0,1,0,0,0,136,82,0,0,0,0,0,0,25,0,0,0,216,234,1,0,1,0,0,0,168,82,0,0,0,0,0,0,25,0,0,0,240,98,2,0,1,0,0,0,168,82,0,0,0,0,0,0,26,0,0,0,144,70,2,0,1,0,0,0,160,82,0,0,0,0,0,0,10,0,0,0,64,55,2,0,1,0,0,0,176,82,0,0,0,0,0,0,11,0,0,0,240,47,2,0,1,0,0,0,176,82,0,0,0,0,0,0,12,0,0,0,240,39,2,0,1,0,0,0,176,82,0,0,0,0,0,0,12,0,0,0,96,33,2,0,1,0,0,0,176,82,0,0,0,0,0,0,12,0,0,0,232,25,2,0,1,0,0,0,176,82,0,0,0,0,0,0,14,0,0,0,0,18,2,0,1,0,0,0,176,82,0,0,0,0,0,0,14,0,0,0,224,7,2,0,1,0,0,0,176,82,0,0,0,0,0,0,13,0,0,0,24,3,2,0,1,0,0,0,176,82,0,0,0,0,0,0,5,0,0,0,144,0,2,0,1,0,0,0,176,82,0,0,0,0,0,0,6,0,0,0,112,254,1,0,1,0,0,0,176,82,0,0,0,0,0,0,7,0,0,0,160,252,1,0,1,0,0,0,176,82,0,0,0,0,0,0,7,0,0,0,136,250,1,0,1,0,0,0,176,82,0,0,0,0,0,0,7,0,0,0,192,247,1,0,1,0,0,0,176,82,0,0,0,0,0,0,9,0,0,0,152,245,1,0,1,0,0,0,176,82,0,0,0,0,0,0,9,0,0,0,184,243,1,0,1,0,0,0,176,82,0,0,0,0,0,0,8,0,0,0,160,240,1,0,1,0,0,0,176,82,0,0,0,0,0,0,0,0,0,0,184,236,1,0,1,0,0,0,128,82,0,0,0,0,0,0,1,0,0,0,40,234,1,0,1,0,0,0,128,82,0,0,0,0,0,0,2,0,0,0,8,232,1,0,1,0,0,0,128,82,0,0,0,0,0,0,2,0,0,0,176,229,1,0,1,0,0,0,128,82,0,0,0,0,0,0,2,0,0,0,8,227,1,0,1,0,0,0,128,82,0,0,0,0,0,0,4,0,0,0,200,224,1,0,1,0,0,0,128,82,0,0,0,0,0,0,4,0,0,0,152,222,1,0,1,0,0,0,128,82,0,0,0,0,0,0,3,0,0,0,160,220,1,0,1,0,0,0,128,82,0,0,0,0,0,0,18,0,0,0,240,95,2,0,1,0,0,0,152,82,0,0,0,0,0,0,27,0,0,0,248,216,1,0,1,0,0,0,144,82,0,0,0,0,0,0,28,0,0,0,192,212,1,0,1,0,0,0,144,82,0,0,0,0,0,0,29,0,0,0,16,210,1,0,1,0,0,0,144,82,0,0,0,0,0,0,29,0,0,0,232,207,1,0,1,0,0,0,144,82,0,0,0,0,0,0,29,0,0,0,248,205,1,0,1,0,0,0,144,82,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,64,206,1,0,0,0,0,0,0,47,0,0,248,46,0,0,1,0,0,0,80,104,2,0,0,0,0,0,120,74,0,0,248,46,0,0,2,0,0,0,104,5,2,0,0,0,0,0,64,15,0,0,248,46,0,0,3,0,0,0,112,235,1,0,0,0,0,0,96,2,0,0,248,46,0,0,4,0,0,0,144,211,1,0,0,0,0,0,88,215,0,0,248,46,0,0,5,0,0,0,24,191,1,0,0,0,0,0,176,34,0,0,248,46,0,0,6,0,0,0,208,170,1,0,0,0,0,0,56,46,0,0,248,46,0,0,7,0,0,0,72,151,1,0,0,0,0,0,160,46,0,0,248,46,0,0,7,0,0,0,144,133,1,0,0,0,0,0,160,46,0,0,248,46,0,0,7,0,0,0,72,117,1,0,0,0,0,0,152,46,0,0,248,46,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,80,235,1,0,0,0,0,0,88,83,0,0,80,83,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,232,243,1,0,78,5,0,0,168,115,2,0,130,4,0,0,128,13,2,0,130,4,0,0,128,238,1,0,84,2,0,0,232,214,1,0,84,2,0,0,240,192,1,0,72,4,0,0,120,173,1,0,72,4,0,0,56,153,1,0,218,1,0,0,8,135,1,0,218,1,0,0,168,118,1,0,32,0,0,0,120,167,2,0,32,0,0,0,88,48,2,0,60,1,0,0,56,138,2,0,60,1,0,0,104,41,2,0,122,3,0,0,0,0,0,0,120,115,1,0,1,0,0,0,0,0,0,0,248,83,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,160,92,2,0,1,0,0,0,0,0,0,0,48,84,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,240,95,2,0,1,0,0,0,0,0,0,0,104,84,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,240,98,2,0,1,0,0,0,0,0,0,0,160,84,0,0,1,0,0,0,248,3,2,0,1,0,0,0,0,0,0,0,160,84,0,0,2,0,0,0,216,234,1,0,1,0,0,0,0,0,0,0,240,85,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,235,1,0,1,0,0,0,0,0,0,0,216,84,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,216,101,2,0,255,255,255,255,0,0,0,0,16,85,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,88,102,2,0,1,0,0,0,0,0,0,0,72,85,0,0,2,0,0,0,160,4,2,0,1,0,0,0,0,0,0,0,128,85,0,0,0,0,0,0,48,235,1,0,1,0,0,0,0,0,0,0,128,85,0,0,3,0,0,0,80,211,1,0,1,0,0,0,0,0,0,0,128,85,0,0,0,0,0,0,216,190,1,0,1,0,0,0,0,0,0,0,72,85,0,0,3,0,0,0,144,170,1,0,1,0,0,0,0,0,0,0,72,85,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,184,102,2,0,1,0,0,0,0,0,0,0,184,85,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,248,190,1,0,1,0,0,0,0,0,0,0,40,86,0,0,0,0,0,0,176,170,1,0,1,0,0,0,0,0,0,0,40,86,0,0,1,0,0,0,24,151,1,0,1,0,0,0,0,0,0,0,96,86,0,0,2,0,0,0,104,133,1,0,1,0,0,0,0,0,0,0,40,86,0,0,3,0,0,0,32,117,1,0,1,0,0,0,0,0,0,0,40,86,0,0,4,0,0,0,168,165,2,0,1,0,0,0,0,0,0,0,40,86,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,120,4,0,0,130,0,0,0,122,5,0,0,198,2,0,0,236,3,0,0,224,4,0,0,26,4,0,0,118,1,0,0,62,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,0,0,0,12,0,0,0,112,3,0,0,0,0,0,0,48,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,96,248,1,0,56,246,1,0,136,218,1,0,0,0,0,0,100,0,0,0,101,0,0,0,102,0,0,0,100,0,0,0,8,241,1,0,160,197,1,0,80,191,1,0,0,0,0,0,0,0,0,0,1,0,0,0,2,0,0,0,255,255,255,255,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,68,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,68,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,65,66,0,0,64,70,2,0,65,73,0,0,16,55,2,0,65,82,0,0,112,47,2,0,65,88,0,0,168,39,2,0,66,32,0,0,64,33,2,0,66,73,0,0,128,25,2,0,67,66,0,0,200,17,2,0,67,79,0,0,144,7,2,0,67,88,0,0,208,2,2,0,72,32,0,0,64,0,2,0,72,66,0,0,64,254,1,0,72,73,0,0,104,252,1,0,72,88,0,0,88,250,1,0,72,98,0,0,136,247,1,0,72,105,0,0,88,245,1,0,72,114,0,0,128,243,1,0,72,120,0,0,104,240,1,0,73,32,0,0,152,236,1,0,75,66,0,0,240,233,1,0,75,73,0,0,200,231,1,0,75,82,0,0,136,229,1,0,75,88,0,0,224,226,1,0,78,66,0,0,160,224,1,0,78,73,0,0,112,222,1,0,78,82,0,0,120,220,1,0,78,88,0,0,232,218,1,0,80,65,0,0,216,216,1,0,80,66,0,0,144,212,1,0,80,73,0,0,240,209,1,0,80,88,0,0,192,207,1,0,82,32,0,0,208,205,1,0,83,32,0,0,8,204,1,0,90,68,0,0,168,202,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,152,149,2,0,200,135,2,0,24,25,2,0,104,94,2,0,16,129,1,0,128,32,2,0,192,6,2,0,104,2,2,0,0,0,0,0,0,0,0,0,126,2,0,0,242,0,0,0,0,0,0,0,0,0,0,0,38,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,152,2,0,0,32,2,0,0,82,0,0,0,164,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,100,3,0,0,2,4,0,0,58,3,0,0,192,2,0,0,144,1,0,0,54,1,0,0,222,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,158,1,0,0,46,4,0,0,0,0,0,0,0,0,240,191,100,1,0,0,150,4,0,0,162,3,0,0,244,2,0,0,0,0,0,0,30,0,0,0,40,221,1,0,198,0,0,0,64,219,1,0,193,0,0,0,48,217,1,0,194,0,0,0,72,213,1,0,192,0,0,0,192,210,1,0,145,3,0,0,120,208,1,0,197,0,0,0,112,206,1,0,195,0,0,0,168,204,1,0,196,0,0,0,232,202,1,0,146,3,0,0,24,201,1,0,199,0,0,0,176,199,1,0,167,3,0,0,248,197,1,0,33,32,0,0,240,194,1,0,148,3,0,0,24,192,1,0,208,0,0,0,56,190,1,0,201,0,0,0,248,187,1,0,202,0,0,0,144,186,1,0,200,0,0,0,232,184,1,0,149,3,0,0,48,183,1,0,151,3,0,0,128,181,1,0,203,0,0,0,184,179,1,0,147,3,0,0,200,177,1,0,205,0,0,0,152,175,1,0,206,0,0,0,216,172,1,0,204,0,0,0,224,169,1,0,153,3,0,0,48,168,1,0,207,0,0,0,24,166,1,0,154,3,0,0,24,164,1,0,155,3,0,0,128,162,1,0,156,3,0,0,144,160,1,0,209,0,0,0,248,158,1,0,157,3,0,0,88,157,1,0,82,1,0,0,32,155,1,0,211,0,0,0,200,152,1,0,212,0,0,0,96,150,1,0,210,0,0,0,24,148,1,0,169,3,0,0,152,146,1,0,159,3,0,0,176,144,1,0,216,0,0,0,104,143,1,0,213,0,0,0,232,141,1,0,214,0,0,0,152,140,1,0,166,3,0,0,48,205,2,0,160,3,0,0,192,136,1,0,51,32,0,0,136,134,1,0,168,3,0,0,152,132,1,0,161,3,0,0,64,131,1,0,96,1,0,0,120,129,1,0,163,3,0,0,48,128,1,0,222,0,0,0,176,126,1,0,164,3,0,0,16,125,1,0,152,3,0,0,104,123,1,0,218,0,0,0,176,121,1,0,219,0,0,0,240,119,1,0,217,0,0,0,32,118,1,0,165,3,0,0,80,116,1,0,220,0,0,0,160,114,1,0,158,3,0,0,112,113,1,0,221,0,0,0,64,112,1,0,120,1,0,0,32,111,1,0,150,3,0,0,56,173,2,0,225,0,0,0,232,171,2,0,226,0,0,0,176,170,2,0,180,0,0,0,232,168,2,0,230,0,0,0,200,166,2,0,224,0,0,0,8,165,2,0,53,33,0,0,80,131,2,0,177,3,0,0,248,161,2,0,38,0,0,0,184,160,2,0,39,34,0,0,104,159,2,0,32,34,0,0,248,157,2,0,229,0,0,0,224,156,2,0,72,34,0,0,120,155,2,0,227,0,0,0,0,154,2,0,228,0,0,0,88,152,2,0,30,32,0,0,112,150,2,0,178,3,0,0,216,148,2,0,166,0,0,0,144,147,2,0,34,32,0,0,96,146,2,0,41,34,0,0,32,145,2,0,231,0,0,0,80,143,2,0,184,0,0,0,48,142,2,0,162,0,0,0,56,141,2,0,199,3,0,0,240,139,2,0,198,2,0,0,240,137,2,0,99,38,0,0,96,136,2,0,69,34,0,0,56,135,2,0,169,0,0,0,240,133,2,0,181,33,0,0,24,132,2,0,42,34,0,0,240,130,2,0,164,0,0,0,192,129,2,0,211,33,0,0,80,128,2,0,32,32,0,0,248,126,2,0,147,33,0,0,64,125,2,0,176,0,0,0,96,123,2,0,180,3,0,0,208,121,2,0,102,38,0,0,104,120,2,0,247,0,0,0,16,119,2,0,233,0,0,0,8,117,2,0,234,0,0,0,160,115,2,0,232,0,0,0,104,114,2,0,5,34,0,0,72,113,2,0,3,32,0,0,16,112,2,0,2,32,0,0,152,7,2,0,181,3,0,0,152,108,2,0,97,34,0,0,40,107,2,0,183,3,0,0,192,105,2,0,240,0,0,0,104,104,2,0,235,0,0,0,120,103,2,0,172,32,0,0,64,102,2,0,3,34,0,0,72,101,2,0,146,1,0,0,80,100,2,0,0,34,0,0,72,99,2,0,189,0,0,0,144,97,2,0,188,0,0,0,112,96,2,0,190,0,0,0,248,94,2,0,68,32,0,0,200,93,2,0,179,3,0,0,120,92,2,0,101,34,0,0,80,91,2,0,62,0,0,0,16,90,2,0,212,33,0,0,0,89,2,0,148,33,0,0,168,87,2,0,101,38,0,0,128,86,2,0,38,32,0,0,184,84,2,0,237,0,0,0,64,83,2,0,238,0,0,0,208,81,2,0,161,0,0,0,112,80,2,0,236,0,0,0,112,79,2,0,17,33,0,0,144,78,2,0,30,34,0,0,144,148,1,0,43,34,0,0,168,76,2,0,185,3,0,0,192,75,2,0,191,0,0,0,216,74,2,0,8,34,0,0,176,73,2,0,239,0,0,0,248,71,2,0,186,3,0,0,176,69,2,0,208,33,0,0,56,65,2,0,187,3,0,0,40,64,2,0,41,35,0,0,48,63,2,0,171,0,0,0,24,62,2,0,144,33,0,0,208,60,2,0,8,35,0,0,136,59,2,0,28,32,0,0,96,58,2,0,100,34,0,0,248,56,2,0,10,35,0,0,208,55,2,0,23,34,0,0,208,54,2,0,202,37,0,0,88,54,2,0,14,32,0,0,240,53,2,0,57,32,0,0,32,53,2,0,24,32,0,0,96,52,2,0,60,0,0,0,192,51,2,0,175,0,0,0,64,51,2,0,20,32,0,0,160,50,2,0,181,0,0,0,136,49,2,0,183,0,0,0,88,48,2,0,18,34,0,0,40,47,2,0,188,3,0,0,64,46,2,0,7,34,0,0,128,45,2,0,160,0,0,0,168,44,2,0,19,32,0,0,240,43,2,0,96,34,0,0,40,43,2,0,11,34,0,0,168,42,2,0,172,0,0,0,16,42,2,0,9,34,0,0,56,41,2,0,132,34,0,0,72,40,2,0,241,0,0,0,104,39,2,0,189,3,0,0,152,38,2,0,243,0,0,0,8,38,2,0,244,0,0,0,144,37,2,0,83,1,0,0,32,37,2,0,242,0,0,0,152,36,2,0,62,32,0,0,16,36,2,0,201,3,0,0,128,35,2,0,191,3,0,0,120,34,2,0,149,34,0,0])
.concat([184,33,2,0,40,34,0,0,232,32,2,0,170,0,0,0,56,32,2,0,186,0,0,0,176,31,2,0,248,0,0,0,32,31,2,0,245,0,0,0,176,30,2,0,151,34,0,0,8,30,2,0,246,0,0,0,144,29,2,0,182,0,0,0,24,29,2,0,2,34,0,0,216,27,2,0,48,32,0,0,136,26,2,0,165,34,0,0,48,25,2,0,198,3,0,0,112,24,2,0,192,3,0,0,200,23,2,0,214,3,0,0,40,23,2,0,177,0,0,0,176,22,2,0,163,0,0,0,176,21,2,0,50,32,0,0,64,21,2,0,15,34,0,0,152,20,2,0,29,34,0,0,96,19,2,0,200,3,0,0,88,18,2,0,34,0,0,0,120,17,2,0,210,33,0,0,160,16,2,0,26,34,0,0,192,15,2,0,42,35,0,0,128,14,2,0,187,0,0,0,136,13,2,0,146,33,0,0,120,12,2,0,9,35,0,0,232,11,2,0,29,32,0,0,40,11,2,0,28,33,0,0,8,10,2,0,174,0,0,0,192,8,2,0,11,35,0,0,56,7,2,0,193,3,0,0,88,6,2,0,15,32,0,0,112,5,2,0,58,32,0,0,232,4,2,0,25,32,0,0,144,4,2,0,26,32,0,0,56,4,2,0,97,1,0,0,24,4,2,0,197,34,0,0,0,4,2,0,167,0,0,0,160,3,2,0,173,0,0,0,96,3,2,0,195,3,0,0,152,2,2,0,194,3,0,0,80,2,2,0,60,34,0,0,40,2,2,0,96,38,0,0,248,1,2,0,130,34,0,0,216,1,2,0,134,34,0,0,176,1,2,0,17,34,0,0,144,1,2,0,131,34,0,0,120,1,2,0,185,0,0,0,40,1,2,0,178,0,0,0,168,0,2,0,179,0,0,0,16,0,2,0,135,34,0,0,192,255,1,0,223,0,0,0,152,255,1,0,196,3,0,0,128,255,1,0,52,34,0,0,104,255,1,0,184,3,0,0,56,255,1,0,209,3,0,0,32,255,1,0,9,32,0,0,8,255,1,0,254,0,0,0,216,156,2,0,220,2,0,0,128,254,1,0,215,0,0,0,32,254,1,0,34,33,0,0,200,253,1,0,209,33,0,0,160,253,1,0,250,0,0,0,136,253,1,0,145,33,0,0,112,253,1,0,251,0,0,0,72,253,1,0,249,0,0,0,48,253,1,0,168,0,0,0,24,253,1,0,210,3,0,0,224,252,1,0,197,3,0,0,184,252,1,0,252,0,0,0,56,252,1,0,24,33,0,0,216,251,1,0,190,3,0,0,176,251,1,0,253,0,0,0,168,251,1,0,165,0,0,0,128,251,1,0,255,0,0,0,88,251,1,0,182,3,0,0,64,251,1,0,13,32,0,0,40,251,1,0,12,32,0,0,184,4,0,0,0,0,0,0,202,1,0,0,0,0,0,0,214,1,0,0,0,0,0,0,114,2,0,0,0,0,0,0,140,3,0,0,0,0,0,0,236,2,0,0,0,0,0,0,78,2,0,0,0,0,0,0,132,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,0,0,0,48,0,0,0,0,0,0,0,210,1,0,0,72,1,0,0,64,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,90,4,0,0,214,0,0,0,0,0,0,0,0,0,0,0,32,1,0,0,238,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,88,64,0,0,0,0,0,0,88,64,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,88,64,0,0,0,0,0,0,88,64,64,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,82,64,0,0,0,0,0,0,82,64,96,0,0,0,0,0,0,0,0,0,0,0,0,0,66,64,0,0,0,0,0,0,66,64,0,0,0,0,0,32,131,64,0,0,0,0,0,192,136,64,0,0,0,0,0,0,82,64,0,0,0,0,0,0,82,64,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,82,64,0,0,0,0,0,0,82,64,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,82,64,0,0,0,0,0,0,82,64,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,88,64,0,0,0,0,0,0,88,64,0,0,12,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,88,64,0,0,0,0,0,0,88,64,2,16,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,128,150,64,0,0,0,0,0,128,150,64,0,0,0,0,0,0,0,0,0,0,0,0,0,0,66,64,0,0,0,0,0,0,66,64,0,0,0,0,0,32,131,64,0,0,0,0,0,192,136,64,0,0,0,0,0,0,82,64,0,0,0,0,0,0,82,64,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,82,64,0,0,0,0,0,0,82,64,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,82,64,0,0,0,0,0,0,82,64,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,48,170,1,0,88,168,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,81,218,27,124,97,50,227,63,12,0,0,0,4,0,0,0,6,0,0,0,2,0,0,0,3,0,0,0,1,0,0,0,9,0,0,0,8,0,0,0,11,0,0,0,12,0,0,0,13,0,0,0,14,0,0,0,15,0,0,0,16,0,0,0,17,0,0,0,18,0,0,0,21,0,0,0,22,0,0,0,23,0,0,0,24,0,0,0,25,0,0,0,26,0,0,0,27,0,0,0,28,0,0,0,31,0,0,0,32,0,0,0,33,0,0,0,34,0,0,0,35,0,0,0,36,0,0,0,37,0,0,0,38,0,0,0,41,0,0,0,42,0,0,0,43,0,0,0,44,0,0,0,45,0,0,0,46,0,0,0,47,0,0,0,48,0,0,0,51,0,0,0,52,0,0,0,53,0,0,0,54,0,0,0,55,0,0,0,56,0,0,0,57,0,0,0,58,0,0,0,61,0,0,0,62,0,0,0,63,0,0,0,64,0,0,0,65,0,0,0,66,0,0,0,67,0,0,0,68,0,0,0,71,0,0,0,72,0,0,0,73,0,0,0,74,0,0,0,75,0,0,0,76,0,0,0,77,0,0,0,78,0,0,0,81,0,0,0,82,0,0,0,83,0,0,0,84,0,0,0,85,0,0,0,86,0,0,0,87,0,0,0,88,0,0,0,8,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,16,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,24,1,2,0,85,93,201,127,201,127,255,0,64,232,1,0,187,45,212,190,174,212,255,0,104,208,1,0,20,119,253,253,192,134,255,0,232,187,1,0,85,93,201,127,201,127,255,0,32,168,1,0,187,45,212,190,174,212,255,0,8,148,1,0,20,119,253,253,192,134,255,0,32,131,1,0,42,102,255,255,255,153,255,0,128,114,1,0,85,93,201,127,201,127,255,0,152,163,2,0,187,45,212,190,174,212,255,0,184,148,2,0,20,119,253,253,192,134,255,0,24,135,2,0,42,102,255,255,255,153,255,0,72,120,2,0,151,173,176,56,108,176,255,0,160,105,2,0,85,93,201,127,201,127,255,0,160,93,2,0,187,45,212,190,174,212,255,0,80,80,2,0,20,119,253,253,192,134,255,0,24,65,2,0,42,102,255,255,255,153,255,0,56,54,2,0,151,173,176,56,108,176,255,0,32,46,2,0,232,252,240,240,2,127,255,0,128,38,2,0,85,93,201,127,201,127,255,0,24,32,2,0,187,45,212,190,174,212,255,0,64,24,2,0,20,119,253,253,192,134,255,0,112,16,2,0,42,102,255,255,255,153,255,0,56,6,2,0,151,173,176,56,108,176,255,0,64,2,2,0,232,252,240,240,2,127,255,0,176,255,1,0,17,224,191,191,91,23,255,0,184,253,1,0,85,93,201,127,201,127,255,0,200,251,1,0,187,45,212,190,174,212,255,0,112,249,1,0,20,119,253,253,192,134,255,0,8,247,1,0,42,102,255,255,255,153,255,0,224,244,1,0,151,173,176,56,108,176,255,0,216,242,1,0,232,252,240,240,2,127,255,0,240,239,1,0,17,224,191,191,91,23,255,0,168,235,1,0,0,0,102,102,102,102,255,0,56,233,1,0,147,25,247,222,235,247,255,0,80,231,1,0,142,75,225,158,202,225,255,0,48,229,1,0,145,188,189,49,130,189,255,0,120,226,1,0,159,16,255,239,243,255,255,0,8,224,1,0,143,46,231,189,215,231,255,0,0,222,1,0,143,127,214,107,174,214,255,0,24,220,1,0,147,208,181,33,113,181,255,0,96,218,1,0,159,16,255,239,243,255,255,0,96,216,1,0,143,46,231,189,215,231,255,0,184,211,1,0,143,127,214,107,174,214,255,0,128,209,1,0,145,188,189,49,130,189,255,0,96,207,1,0,149,241,156,8,81,156,255,0,80,205,1,0,159,16,255,239,243,255,255,0,168,203,1,0,148,43,239,198,219,239,255,0,48,202,1,0,142,75,225,158,202,225,255,0,88,200,1,0,143,127,214,107,174,214,255,0,168,198,1,0,145,188,189,49,130,189,255,0,232,196,1,0,149,241,156,8,81,156,255,0,208,193,1,0,159,16,255,239,243,255,255,0,64,191,1,0,148,43,239,198,219,239,255,0,32,189,1,0,142,75,225,158,202,225,255,0,200,188,1,0,143,127,214,107,174,214,255,0,0,187,1,0,144,169,198,66,146,198,255,0,104,185,1,0,147,208,181,33,113,181,255,0,184,183,1,0,151,241,148,8,69,148,255,0,120,180,1,0,148,8,255,247,251,255,255,0,112,178,1,0,147,25,247,222,235,247,255,0,48,178,1,0,148,43,239,198,219,239,255,0,96,176,1,0,142,75,225,158,202,225,255,0,104,173,1,0,143,127,214,107,174,214,255,0,120,170,1,0,144,169,198,66,146,198,255,0,152,168,1,0,147,208,181,33,113,181,255,0,136,166,1,0,151,241,148,8,69,148,255,0,200,164,1,0,148,8,255,247,251,255,255,0,32,163,1,0,147,25,247,222,235,247,255,0,40,161,1,0,148,43,239,198,219,239,255,0,112,159,1,0,142,75,225,158,202,225,255,0,200,157,1,0,143,127,214,107,174,214,255,0,240,155,1,0,144,169,198,66,146,198,255,0,40,153,1,0,147,208,181,33,113,181,255,0,216,150,1,0,149,241,156,8,81,156,255,0,128,148,1,0,152,235,107,8,48,107,255,0,240,146,1,0,23,239,84,84,48,5,255,0,24,145,1,0,119,255,60,0,60,48,255,0,192,143,1,0,23,236,140,140,81,10,255,0,72,142,1,0,24,194,191,191,129,45,255,0,248,140,1,0,29,112,223,223,194,125,255,0,144,139,1,0,30,52,246,246,232,195,255,0,232,137,1,0,121,38,234,199,234,229,255,0,32,135,1,0,120,95,205,128,205,193,255,0,240,131,1,0,124,165,151,53,151,143,255,0,32,130,1,0,124,252,102,1,102,94,255,0,224,129,1,0,23,239,84,84,48,5,255,0,136,128,1,0,124,252,102,1,102,94,255,0,24,127,1,0,119,255,60,0,60,48,255,0,128,125,1,0,23,236,140,140,81,10,255,0,216,123,1,0,24,194,191,191,129,45,255,0,56,122,1,0,29,112,223,223,194,125,255,0,184,120,1,0,30,52,246,246,232,195,255,0,192,118,1,0,0,0,245,245,245,245,255,0,240,116,1,0,121,38,234,199,234,229,255,0,24,115,1,0,120,95,205,128,205,193,255,0,216,113,1,0,124,165,151,53,151,143,255,0,152,112,1,0,28,135,216,216,179,101,255,0,136,111,1,0,0,0,245,245,245,245,255,0,208,173,2,0,123,127,180,90,180,172,255,0,88,172,2,0,21,215,166,166,97,26,255,0,72,171,2,0,29,112,223,223,194,125,255,0,208,169,2,0,120,95,205,128,205,193,255,0,144,167,2,0,121,253,133,1,133,113,255,0,120,165,2,0,21,215,166,166,97,26,255,0,160,162,2,0,29,112,223,223,194,125,255,0,96,162,2,0,0,0,245,245,245,245,255,0,0,161,2,0,120,95,205,128,205,193,255,0,176,159,2,0,121,253,133,1,133,113,255,0,64,158,2,0,23,236,140,140,81,10,255,0,40,157,2,0,28,135,216,216,179,101,255,0,248,155,2,0,30,52,246,246,232,195,255,0,120,154,2,0,121,38,234,199,234,229,255,0,184,152,2,0,123,127,180,90,180,172,255,0,32,151,2,0,124,252,102,1,102,94,255,0,96,149,2,0,23,236,140,140,81,10,255,0,8,148,2,0,28,135,216,216,179,101,255,0,184,146,2,0,30,52,246,246,232,195,255,0,152,145,2,0,0,0,245,245,245,245,255,0,208,143,2,0,121,38,234,199,234,229,255,0,144,142,2,0,123,127,180,90,180,172,255,0,136,141,2,0,124,252,102,1,102,94,255,0,88,140,2,0,23,236,140,140,81,10,255,0,104,138,2,0,24,194,191,191,129,45,255,0,184,136,2,0,29,112,223,223,194,125,255,0,144,135,2,0,30,52,246,246,232,195,255,0,248,132,2,0,121,38,234,199,234,229,255,0,120,132,2,0,120,95,205,128,205,193,255,0,64,131,2,0,124,165,151,53,151,143,255,0,40,130,2,0,124,252,102,1,102,94,255,0,216,128,2,0,23,236,140,140,81,10,255,0,112,127,2,0,24,194,191,191,129,45,255,0,152,125,2,0,29,112,223,223,194,125,255,0,200,123,2,0,30,52,246,246,232,195,255,0,24,122,2,0,0,0,245,245,245,245,255,0,176,120,2,0,121,38,234,199,234,229,255,0,96,119,2,0,120,95,205,128,205,193,255,0,248,117,2,0,124,165,151,53,151,143,255,0,48,116,2,0,124,252,102,1,102,94,255,0,200,114,2,0,135,20,249,229,245,249,255,0,144,113,2,0,117,74,216,153,216,201,255,0,88,112,2,0,103,185,162,44,162,95,255,0,8,111,2,0,136,14,251,237,248,251,255,0,16,109,2,0,127,54,226,178,226,226,255,0,112,107,2,0,113,120,194,102,194,164,255,0,48,106,2,0,98,190,139,35,139,69,255,0,192,104,2,0,136,14,251,237,248,251,255,0,192,103,2,0,127,54,226,178,226,226,255,0,152,102,2,0,113,120,194,102,194,164,255,0,144,101,2,0,103,185,162,44,162,95,255,0,152,100,2,0,102,255,109,0,109,44,255,0,120,98,2,0,136,14,251,237,248,251,255,0,48,98,2,0,119,34,236,204,236,230,255,0,208,96,2,0,117,74,216,153,216,201,255,0,104,95,2,0,113,120,194,102,194,164,255,0,24,93,2,0,103,185,162,44,162,95,255,0,208,92,2,0,102,255,109,0,109,44,255,0,152,91,2,0,136,14,251,237,248,251,255,0,112,90,2,0,119,34,236,204,236,230,255,0,104,89,2,0,117,74,216,153,216,201,255,0,16,88,2,0,113,120,194,102,194,164,255,0,168,85,2,0,105,159,174,65,174,118,255,0,104,85,2,0,98,190,139,35,139,69,255,0,184,83,2,0,102,255,88,0,88,36,255,0,80,82,2,0,134,6,253,247,252,253,255,0,192,80,2,0,135,20,249,229,245,249,255,0,184,79,2,0,119,34,236,204,236,230,255,0,216,78,2,0,117,74,216,153,216,201,255,0,208,77,2,0,113,120,194,102,194,164,255,0,240,76,2,0,105,159,174,65,174,118,255,0,8,76,2,0,98,190,139,35,139,69,255,0,32,75,2,0,102,255,88,0,88,36,255,0,32,74,2,0,134,6,253,247,252,253,255,0,112,72,2,0,135,20,249,229,245,249,255,0,48,70,2,0,119,34,236,204,236,230,255,0,152,65,2,0,117,74,216,153,216,201,255,0,112,64,2,0,113,120,194,102,194,164,255,0,160,62,2,0,105,159,174,65,174,118,255,0,96,62,2,0,98,190,139,35,139,69,255,0,24,61,2,0,102,255,109,0,109,44,255,0,208,59,2,0,101,255,68,0,68,27,255,0,168,58,2,0,144,20,244,224,236,244,255,0,160,57,2,0,148,70,218,158,188,218,255,0,72,56,2,0,196,123,167,136,86,167,255,0,0,55,2,0,136,14,251,237,248,251,255,0,104,54,2,0,146,53,227,179,205,227,255,0,136,53,2,0,162,74,198,140,150,198,255,0,88,53,2,0,202,149,157,136,65,157,255,0,152,52,2,0,136,14,251,237,248,251,255,0,8,52,2,0,146,53,227,179,205,227,255,0,120,51,2,0,162,74,198,140,150,198,255,0,216,50,2,0,196,123,167,136,86,167,255,0,248,49,2,0,214,225,129,129,15,124,255,0,248,47,2,0,136,14,251,237,248,251,255,0,96,47,2,0,148,43,230,191,211,230,255,0,136,46,2,0,148,70,218,158,188,218,255,0,136,45,2,0,162,74,198,140,150,198,255,0,224,44,2,0,196,123,167,136,86,167,255,0,112,43,2,0,214,225,129,129,15,124,255,0,72,43,2,0,136,14,251,237,248,251,255,0,208,42,2,0,148,43,230,191,211,230,255,0,48,42,2,0,148,70,218,158,188,218,255,0,128,41,2,0,162,74,198,140,150,198,255,0,200,40,2,0,190,100,177,140,107,177,255,0,152,39,2,0,202,149,157,136,65,157,255,0,192,38,2,0,213,252,110,110,1,107,255,0,48,38,2,0,134,6,253,247,252,253,255,0,176,37,2,0,144,20,244,224,236,244,255,0,56,37,2,0,148,43,230,191,211,230,255,0,192,36,2,0,148,70,218,158,188,218,255,0,56,36,2,0,162,74,198,140,150,198,255,0,168,35,2,0,190,100,177,140,107,177,255,0,48,35,2,0,202,149,157,136,65,157,255,0,248,33,2,0,213,252,110,110,1,107,255,0,48,33,2,0,134,6,253,247,252,253,255,0,96,32,2,0,144,20,244,224,236,244,255,0,216,31,2,0,148,43,230,191,211,230,255,0,96,31,2,0,148,70,218,158,188,218,255,0,208,30,2,0,162,74,198,140,150,198,255,0,64,30,2,0,190,100,177,140,107,177,255,0,176,29,2,0,202,149,157,136,65,157,255,0,64,29,2,0,214,225,129,129,15,124,255,0,144,27,2,0,213,255,77,77,0,75,255,0,64,27,2,0,114,211,158,27,158,119,255,0,112,25,2,0,18,252,217,217,95,2,255,0,160,24,2,0,173,95,179,117,112,179,255,0,0,24,2,0,114,211,158,27,158,119,255,0,64,23,2,0,18,252,217,217,95,2,255,0,200,22,2,0,173,95,179,117,112,179,255,0,24,22,2,0,233,209,231,231,41,138,255,0,128,21,2,0,114,211,158,27,158,119,255,0,192,20,2,0,18,252,217,217,95,2,255,0,56,20,2,0,173,95,179,117,112,179,255,0,224,18,2,0,233,209,231,231,41,138,255,0,184,17,2,0,62,208,166,102,166,30,255,0,200,16,2,0,114,211,158,27,158,119,255,0,248,15,2,0,18,252,217,217,95,2,255,0,40,15,2,0,173,95,179,117,112,179,255,0,248,13,2,0,233,209,231,231,41,138,255,0,216,12,2,0,62,208,166,102,166,30,255,0,32,12,2,0,31,252,230,230,171,2,255,0,80,11,2,0,114,211,158,27,158,119,255,0,168,10,2,0,18,252,217,217,95,2,255,0,32,9,2,0,173,95,179,117,112,179,255,0,128,7,2,0,233,209,231,231,41,138,255,0,152,6,2,0,62,208,166,102,166,30,255,0,192,5,2,0,31,252,230,230,171,2,255,0,8,5,2,0,27,210,166,166,118,29,255,0,192,4,2,0,114,211,158,27,158,119,255,0,64,4,2,0,18,252,217,217,95,2,255,0,32,4,2,0,173,95,179,117,112,179,255,0,8,4,2,0,233,209,231,231,41,138,255,0,216,3,2,0,62,208,166,102,166,30,255,0,112,3,2,0,31,252,230,230,171,2,255,0,192,2,2,0,27,210,166,166,118,29,255,0,88,2,2,0,0,0,102,102,102,102,255,0,16,2,2,0,76,25,243,224,243,219,255,0,0,2,2,0,95,61,221,168,221,181,255,0,232,1,2,0,140,170,202,67,162,202,255,0,184,1,2,0,65,17,249,240,249,232,255,0,160,1,2,0,87,46,228,186,228,188,255,0,128,1,2,0,123,101,204,123,204,196,255,0,80,1,2,0,141,197,190,43,140,190,255,0,232,0,2,0,65,17,249,240,249,232,255,0,48,0,2,0,87,46,228,186,228,188,255,0,200,255,1,0,123,101,204,123,204,196,255,0,160,255,1,0,140,170,202,67,162,202,255,0,136,255,1,0,145,243,172,8,104,172,255,0,112,255,1,0,65,17,249,240,249,232,255,0,72,255,1,0,77,41,235,204,235,197,255,0,40,255,1,0,95,61,221,168,221,181,255,0,16,255,1,0,123,101,204,123,204,196,255,0,216,254,1,0,140,170,202,67,162,202,255,0,176,254,1,0,145,243,172,8,104,172,255,0,248,253,1,0,65,17,249,240,249,232,255,0,208,253,1,0,77,41,235,204,235,197,255,0,168,253,1,0,95,61,221,168,221,181,255,0,144,253,1,0,123,101,204,123,204,196,255,0,120,253,1,0,137,160,211,78,179,211,255,0,80,253,1,0,141,197,190,43,140,190,255,0,56,253,1,0,147,242,158,8,88,158,255,0,32,253,1,0,60,12,252,247,252,240,255,0,240,252,1,0,76,25,243,224,243,219,255,0,208,252,1,0,77,41,235,204,235,197,255,0,88,252,1,0,95,61,221,168,221,181,255,0,224,251,1,0,123,101,204,123,204,196,255,0,184,251,1,0,137,160,211,78,179,211,255,0,152,251,1,0,141,197,190,43,140,190,255,0,136,251,1,0,147,242,158,8,88,158,255,0,96,251,1,0,60,12,252,247,252,240,255,0,72,251,1,0,76,25,243,224,243,219,255,0,48,251,1,0,77,41,235,204,235,197,255,0,24,251,1,0,95,61,221,168,221,181,255,0,232,250,1,0,123,101,204,123,204,196,255,0,72,250,1,0,137,160,211,78,179,211,255,0,168,249,1,0,141,197,190,43,140,190,255,0,96,249,1,0,145,243,172,8,104,172,255,0,40,249,1,0,150,239,129,8,64,129,255,0,240,248,1,0,74,21,245,229,245,224,255,0,200,248,1,0,80,72,217,161,217,155,255,0,184,248,1,0,98,178,163,49,163,84,255,0,136,248,1,0,73,15,248,237,248,233,255,0,80,248,1,0,78,54,228,186,228,179,255,0,40,248,1,0,86,104,196,116,196,118,255,0,120,247,1,0,98,190,139,35,139,69,255,0,32,247,1,0,73,15,248,237,248,233,255,0,240,246,1,0,78,54,228,186,228,179,255,0,224,246,1,0,86,104,196,116,196,118,255,0,208,246,1,0,98,178,163,49,163,84,255,0,176,246,1,0,102,255,109,0,109,44,255,0,160,246,1,0,73,15,248,237,248,233,255,0,144,246,1,0,77,44,233,199,233,192,255,0,40,246,1,0,80,72,217,161,217,155,255,0,248,245,1,0,86,104,196,116,196,118,255,0,72,245,1,0,98,178,163,49,163,84,255,0,240,244,1,0,102,255,109,0,109,44,255,0,208,244,1,0,73,15,248,237,248,233,255,0,192,244,1,0,77,44,233,199,233,192,255,0,152,244,1,0,80,72,217,161,217,155,255,0,120,244,1,0,86,104,196,116,196,118,255,0,104,244,1,0,96,158,171,65,171,93,255,0,88,244,1,0,98,190,139,35,139,69,255,0,248,243,1,0,108,255,90,0,90,50,255,0,208,243,1,0,72,7,252,247,252,245,255,0,112,243,1,0,74,21,245,229,245,224,255,0,232,242,1,0,77,44,233,199,233,192,255,0,200,242,1,0,80,72,217,161,217,155,255,0,136,242,1,0,86,104,196,116,196,118,255,0,120,242,1,0,96,158,171,65,171,93,255,0,96,242,1,0,98,190,139,35,139,69,255,0,32,242,1,0,108,255,90,0,90,50,255,0,144,241,1,0,72,7,252,247,252,245,255,0,248,240,1,0,74,21,245,229,245,224,255,0,176,240,1,0,77,44,233,199,233,192,255,0,88,240,1,0,80,72,217,161,217,155,255,0,0,240,1,0,86,104,196,116,196,118,255,0,184,239,1,0,96,158,171,65,171,93,255,0,24,239,1,0,98,190,139,35,139,69,255,0,176,238,1,0,102,255,109,0,109,44,255,0,96,238,1,0,101,255,68,0,68,27,255,0,0,238,1,0,0,0,240,240,240,240,255,0,232,237,1,0,0,0,189,189,189,189,255,0,184,237,1,0,0,0,99,99,99,99,255,0,96,237,1,0,0,0,247,247,247,247,255,0,136,236,1,0,0,0,204,204,204,204,255,0,224,235,1,0,0,0,150,150,150,150,255,0,144,235,1,0,0,0,82,82,82,82,255,0,88,235,1,0,0,0,247,247,247,247,255,0,64,235,1,0,0,0,204,204,204,204,255,0,8,235,1,0,0,0,150,150,150,150,255,0,240,234,1,0,0,0,99,99,99,99,255,0,224,234,1,0,0,0,37,37,37,37,255,0,184,234,1,0,0,0,247,247,247,247,255,0,104,234,1,0,0,0,217,217,217,217,255,0,0,234,1,0,0,0,189,189,189,189,255,0,160,233,1,0,0,0,150,150,150,150,255,0,40,233,1,0,0,0,99,99,99,99,255,0,16,233,1,0,0,0,37,37,37,37,255,0,0,233,1,0,0,0,247,247,247,247,255,0,216,232,1,0,0,0,217,217,217,217,255,0,200,232,1,0,0,0,189,189,189,189,255,0,184,232,1,0,0,0,150,150,150,150,255,0,88,232,1,0,0,0,115,115,115,115,255,0,48,232,1,0,0,0,82,82,82,82,255,0,224,231,1,0,0,0,37,37,37,37,255,0,112,231,1,0,0,0,255,255,255,255,255,0,32,231,1,0,0,0,240,240,240,240,255,0,16,231,1,0,0,0,217,217,217,217,255,0,216,230,1,0,0,0,189,189,189,189,255,0,192,230,1,0,0,0,150,150,150,150,255,0,176,230,1,0,0,0,115,115,115,115,255,0,160,230,1,0,0,0,82,82,82,82,255,0,64,230,1,0,0,0,37,37,37,37,255,0,16,230,1,0,0,0,255,255,255,255,255,0,152,229,1,0,0,0,240,240,240,240,255,0,32,229,1,0,0,0,217,217,217,217,255,0,232,228,1,0,0,0,189,189,189,189,255,0,216,228,1,0,0,0,150,150,150,150,255,0,200,228,1,0,0,0,115,115,115,115,255,0,80,228,1,0,0,0,82,82,82,82,255,0,56,228,1,0,0,0,37,37,37,37,255,0,40,228,1,0,0,0,0,0,0,0,255,0,184,227,1,0,21,48,254,254,230,206,255,0,152,227,1,0,19,147,253,253,174,107,255,0,248,226,1,0,14,240,230,230,85,13,255,0,136,226,1,0,19,32,254,254,237,222,255,0,104,226,1,0,20,120,253,253,190,133,255,0,56,226,1,0,17,194,253,253,141,60,255,0,40,226,1,0,13,253,217,217,71,1,255,0,8,226,1,0,19,32,254,254,237,222,255,0,200,225,1,0,20,120,253,253,190,133,255,0,152,225,1,0,17,194,253,253,141,60,255,0,120,225,1,0,14,240,230,230,85,13,255,0,32,225,1,0,13,250,166,166,54,3,255,0,184,224,1,0,19,32,254,254,237,222,255,0,64,224,1,0,21,91,253,253,208,162,255,0,248,223,1,0,19,147,253,253,174,107,255,0,232,223,1,0,17,194,253,253,141,60,255,0,216,223,1,0,14,240,230,230,85,13,255,0,152,223,1,0,13,250,166,166,54,3,255,0,40,223,1,0,19,32,254,254,237,222,255,0,24,223,1,0,21,91,253,253,208,162,255,0,224,222,1,0,19,147,253,253,174,107,255,0,184,222,1,0,17,194,253,253,141,60,255,0,136,222,1,0,16,234,241,241,105,19,255,0,16,222,1,0,13,253,217,217,72,1,255,0,240,221,1,0,12,247,140,140,45,4,255,0,208,221,1,0,21,20,255,255,245,235,255,0,192,221,1,0,21,48,254,254,230,206,255,0,160,221,1,0,21,91,253,253,208,162,255,0,144,221,1,0,19,147,253,253,174,107,255,0,128,221,1,0,17,194,253,253,141,60,255,0,64,221,1,0,16,234,241,241,105,19,255,0,232,220,1,0,13,253,217,217,72,1,255,0,144,220,1,0,12,247,140,140,45,4,255,0,40,220,1,0,21,20,255,255,245,235,255,0,8,220,1,0,21,48,254,254,230,206,255,0,248,219,1,0,21,91,253,253,208,162,255,0,232,219,1,0,19,147,253,253,174,107,255,0,200,219,1,0,17,194,253,253,141,60,255,0,184,219,1,0,16,234,241,241,105,19,255,0,168,219,1,0,13,253,217,217,72,1,255,0,88,219,1,0,13,250,166,166,54,3,255,0,48,219,1,0,12,246,127,127,39,4,255,0,8,219,1,0,25,54,254,254,232,200,255,0,120,218,1,0,19,121,253,253,187,132,255,0,80,218,1,0,5,197,227,227,74,51,255,0,56,218,1,0,26,37,254,254,240,217,255,0,40,218,1,0,24,115,253,253,204,138,255,0,0,218,1,0,13,164,252,252,141,89,255,0,240,217,1,0,3,218,215,215,48,31,255,0,144,217,1,0,26,37,254,254,240,217,255,0,80,217,1,0,24,115,253,253,204,138,255,0,32,217,1,0,13,164,252,252,141,89,255,0,232,216,1,0,5,197,227,227,74,51,255,0,112,216,1,0,0,255,179,179,0,0,255,0,48,216,1,0,26,37,254,254,240,217,255,0,16,216,1,0,24,95,253,253,212,158,255,0,184,214,1,0,19,121,253,253,187,132,255,0,152,214,1,0,13,164,252,252,141,89,255,0,16,214,1,0,5,197,227,227,74,51,255,0,0,214,1,0,0,255,179,179,0,0,255,0,136,213,1,0,26,37,254,254,240,217,255,0,24,213,1,0,24,95,253,253,212,158,255,0,160,212,1,0,19,121,253,253,187,132,255,0,248,211,1,0,13,164,252,252,141,89,255,0,152,211,1,0,7,178,239,239,101,72,255,0,120,211,1,0,3,218,215,215,48,31,255,0,96,211,1,0,0,255,153,153,0,0,255,0,40,211,1,0,24,18,255,255,247,236,255,0,24,211,1,0,25,54,254,254,232,200,255,0,8,211,1,0,24,95,253,253,212,158,255,0,208,210,1,0,19,121,253,253,187,132,255,0,152,210,1,0,13,164,252,252,141,89,255,0,0,210,1,0,7,178,239,239,101,72,255,0,176,209,1,0,3,218,215,215,48,31,255,0,112,209,1,0,0,255,153,153,0,0,255,0,72,209,1,0,24,18,255,255,247,236,255,0,0,209,1,0,25,54,254,254,232,200,255,0,224,208,1,0,24,95,253,253,212,158,255,0,200,208,1,0,19,121,253,253,187,132,255,0,184,208,1,0,13,164,252,252,141,89,255,0,144,208,1,0,7,178,239,239,101,72,255,0,88,208,1,0,3,218,215,215,48,31,255,0,216,207,1,0,0,255,179,179,0,0,255,0,112,207,1,0,0,255,127,127,0,0,255,0,80,207,1,0,142,68,227,166,206,227,255,0,64,207,1,0,190,153,154,106,61,154,255,0,8,207,1,0,144,211,180,31,120,180,255,0,224,206,1,0,65,97,223,178,223,138,255,0,208,206,1,0,82,184,160,51,160,44,255,0,192,206,1,0,0,99,251,251,154,153,255,0,128,206,1,0,254,225,227,227,26,28,255,0,80,206,1,0,23,143,253,253,191,111,255,0,224,205,1,0,21,255,255,255,127,0,255,0,96,205,1,0,198,42,214,202,178,214,255,0,64,205,1,0,142,68,227,166,206,227,255,0,48,205,1,0,190,153,154,106,61,154,255,0,32,205,1,0,42,102,255,255,255,153,255,0,8,205,1,0,144,211,180,31,120,180,255,0,248,204,1,0,65,97,223,178,223,138,255,0,232,204,1,0,82,184,160,51,160,44,255,0,192,204,1,0,0,99,251,251,154,153,255,0,144,204,1,0,254,225,227,227,26,28,255,0,16,204,1,0,23,143,253,253,191,111,255,0,184,203,1,0,21,255,255,255,127,0,255,0,152,203,1,0,198,42,214,202,178,214,255,0,128,203,1,0,142,68,227,166,206,227,255,0,112,203,1,0,190,153,154,106,61,154,255,0,80,203,1,0,42,102,255,255,255,153,255,0,64,203,1,0,15,197,177,177,89,40,255,0,48,203,1,0,144,211,180,31,120,180,255,0,248,202,1,0,65,97,223,178,223,138,255,0,216,202,1,0,82,184,160,51,160,44,255,0,184,202,1,0,0,99,251,251,154,153,255,0,72,202,1,0,254,225,227,227,26,28,255,0,32,202,1,0,23,143,253,253,191,111,255,0,16,202,1,0,21,255,255,255,127,0,255,0,0,202,1,0,198,42,214,202,178,214,255,0,208,201,1,0,142,68,227,166,206,227,255,0,192,201,1,0,144,211,180,31,120,180,255,0,176,201,1,0,65,97,223,178,223,138,255,0,48,201,1,0,142,68,227,166,206,227,255,0,8,201,1,0,144,211,180,31,120,180,255,0,128,200,1,0,65,97,223,178,223,138,255,0,104,200,1,0,82,184,160,51,160,44,255,0,72,200,1,0,142,68,227,166,206,227,255,0,56,200,1,0,144,211,180,31,120,180,255,0,40,200,1,0,65,97,223,178,223,138,255,0,8,200,1,0,82,184,160,51])
.concat([160,44,255,0,248,199,1,0,0,99,251,251,154,153,255,0,232,199,1,0,142,68,227,166,206,227,255,0,200,199,1,0,144,211,180,31,120,180,255,0,160,199,1,0,65,97,223,178,223,138,255,0,112,199,1,0,82,184,160,51,160,44,255,0,200,198,1,0,0,99,251,251,154,153,255,0,152,198,1,0,254,225,227,227,26,28,255,0,136,198,1,0,142,68,227,166,206,227,255,0,120,198,1,0,144,211,180,31,120,180,255,0,88,198,1,0,65,97,223,178,223,138,255,0,72,198,1,0,82,184,160,51,160,44,255,0,56,198,1,0,0,99,251,251,154,153,255,0,8,198,1,0,254,225,227,227,26,28,255,0,232,197,1,0,23,143,253,253,191,111,255,0,136,197,1,0,142,68,227,166,206,227,255,0,24,197,1,0,144,211,180,31,120,180,255,0,200,196,1,0,65,97,223,178,223,138,255,0,160,196,1,0,82,184,160,51,160,44,255,0,136,196,1,0,0,99,251,251,154,153,255,0,224,195,1,0,254,225,227,227,26,28,255,0,208,195,1,0,23,143,253,253,191,111,255,0,64,195,1,0,21,255,255,255,127,0,255,0,32,195,1,0,142,68,227,166,206,227,255,0,224,194,1,0,144,211,180,31,120,180,255,0,192,194,1,0,65,97,223,178,223,138,255,0,224,193,1,0,82,184,160,51,160,44,255,0,168,193,1,0,0,99,251,251,154,153,255,0,152,193,1,0,254,225,227,227,26,28,255,0,72,193,1,0,23,143,253,253,191,111,255,0,216,192,1,0,21,255,255,255,127,0,255,0,192,192,1,0,198,42,214,202,178,214,255,0,120,192,1,0,3,78,251,251,180,174,255,0,56,192,1,0,146,53,227,179,205,227,255,0,8,192,1,0,77,41,235,204,235,197,255,0,216,191,1,0,3,78,251,251,180,174,255,0,112,191,1,0,146,53,227,179,205,227,255,0,48,191,1,0,77,41,235,204,235,197,255,0,0,191,1,0,202,27,228,222,203,228,255,0,232,190,1,0,3,78,251,251,180,174,255,0,200,190,1,0,146,53,227,179,205,227,255,0,168,190,1,0,77,41,235,204,235,197,255,0,136,190,1,0,202,27,228,222,203,228,255,0,80,190,1,0,24,88,254,254,217,166,255,0,232,189,1,0,3,78,251,251,180,174,255,0,192,189,1,0,146,53,227,179,205,227,255,0,48,189,1,0,77,41,235,204,235,197,255,0,16,189,1,0,202,27,228,222,203,228,255,0,240,188,1,0,24,88,254,254,217,166,255,0,224,188,1,0,42,50,255,255,255,204,255,0,184,188,1,0,3,78,251,251,180,174,255,0,168,188,1,0,146,53,227,179,205,227,255,0,152,188,1,0,77,41,235,204,235,197,255,0,16,188,1,0,202,27,228,222,203,228,255,0,216,187,1,0,24,88,254,254,217,166,255,0,176,187,1,0,42,50,255,255,255,204,255,0,64,187,1,0,28,44,229,229,216,189,255,0,48,187,1,0,3,78,251,251,180,174,255,0,32,187,1,0,146,53,227,179,205,227,255,0,16,187,1,0,77,41,235,204,235,197,255,0,240,186,1,0,202,27,228,222,203,228,255,0,224,186,1,0,24,88,254,254,217,166,255,0,200,186,1,0,42,50,255,255,255,204,255,0,168,186,1,0,28,44,229,229,216,189,255,0,128,186,1,0,233,35,253,253,218,236,255,0,96,186,1,0,3,78,251,251,180,174,255,0,200,185,1,0,146,53,227,179,205,227,255,0,152,185,1,0,77,41,235,204,235,197,255,0,136,185,1,0,202,27,228,222,203,228,255,0,120,185,1,0,24,88,254,254,217,166,255,0,88,185,1,0,42,50,255,255,255,204,255,0,72,185,1,0,28,44,229,229,216,189,255,0,56,185,1,0,233,35,253,253,218,236,255,0,0,185,1,0,0,0,242,242,242,242,255,0,216,184,1,0,108,53,226,179,226,205,255,0,152,184,1,0,17,81,253,253,205,172,255,0,32,184,1,0,155,31,232,203,213,232,255,0,232,183,1,0,108,53,226,179,226,205,255,0,216,183,1,0,17,81,253,253,205,172,255,0,200,183,1,0,155,31,232,203,213,232,255,0,168,183,1,0,228,43,244,244,202,228,255,0,144,183,1,0,108,53,226,179,226,205,255,0,104,183,1,0,17,81,253,253,205,172,255,0,72,183,1,0,155,31,232,203,213,232,255,0,32,183,1,0,228,43,244,244,202,228,255,0,240,182,1,0,56,45,245,230,245,201,255,0,120,182,1,0,108,53,226,179,226,205,255,0,80,182,1,0,17,81,253,253,205,172,255,0,64,182,1,0,155,31,232,203,213,232,255,0,48,182,1,0,228,43,244,244,202,228,255,0,232,181,1,0,56,45,245,230,245,201,255,0,216,181,1,0,35,81,255,255,242,174,255,0,200,181,1,0,108,53,226,179,226,205,255,0,152,181,1,0,17,81,253,253,205,172,255,0,112,181,1,0,155,31,232,203,213,232,255,0,32,181,1,0,228,43,244,244,202,228,255,0,160,180,1,0,56,45,245,230,245,201,255,0,104,180,1,0,35,81,255,255,242,174,255,0,88,180,1,0,25,39,241,241,226,204,255,0,72,180,1,0,108,53,226,179,226,205,255,0,24,180,1,0,17,81,253,253,205,172,255,0,8,180,1,0,155,31,232,203,213,232,255,0,248,179,1,0,228,43,244,244,202,228,255,0,208,179,1,0,56,45,245,230,245,201,255,0,168,179,1,0,35,81,255,255,242,174,255,0,128,179,1,0,25,39,241,241,226,204,255,0,128,178,1,0,0,0,204,204,204,204,255,0,96,178,1,0,230,253,142,142,1,82,255,0,80,178,1,0,77,191,100,39,100,25,255,0,64,178,1,0,230,220,197,197,27,125,255,0,32,178,1,0,232,118,222,222,119,174,255,0,16,178,1,0,229,62,241,241,182,218,255,0,0,178,1,0,233,29,253,253,224,239,255,0,224,177,1,0,59,38,245,230,245,208,255,0,184,177,1,0,61,103,225,184,225,134,255,0,120,177,1,0,63,166,188,127,188,65,255,0,200,176,1,0,68,197,146,77,146,33,255,0,160,176,1,0,230,253,142,142,1,82,255,0,144,176,1,0,68,197,146,77,146,33,255,0,120,176,1,0,77,191,100,39,100,25,255,0,80,176,1,0,230,220,197,197,27,125,255,0,8,176,1,0,232,118,222,222,119,174,255,0,184,175,1,0,229,62,241,241,182,218,255,0,168,175,1,0,233,29,253,253,224,239,255,0,128,175,1,0,0,0,247,247,247,247,255,0,88,175,1,0,59,38,245,230,245,208,255,0,208,174,1,0,61,103,225,184,225,134,255,0,144,174,1,0,63,166,188,127,188,65,255,0,128,174,1,0,231,76,233,233,163,201,255,0,48,174,1,0,0,0,247,247,247,247,255,0,88,173,1,0,63,129,215,161,215,106,255,0,16,173,1,0,228,220,208,208,28,139,255,0,0,173,1,0,229,62,241,241,182,218,255,0,240,172,1,0,61,103,225,184,225,134,255,0,200,172,1,0,72,198,172,77,172,38,255,0,112,172,1,0,228,220,208,208,28,139,255,0,112,171,1,0,229,62,241,241,182,218,255,0,48,171,1,0,0,0,247,247,247,247,255,0,184,170,1,0,61,103,225,184,225,134,255,0,160,170,1,0,72,198,172,77,172,38,255,0,104,170,1,0,230,220,197,197,27,125,255,0,88,170,1,0,231,76,233,233,163,201,255,0,72,170,1,0,233,29,253,253,224,239,255,0,248,169,1,0,59,38,245,230,245,208,255,0,168,169,1,0,63,129,215,161,215,106,255,0,112,169,1,0,68,197,146,77,146,33,255,0,0,169,1,0,230,220,197,197,27,125,255,0,240,168,1,0,231,76,233,233,163,201,255,0,200,168,1,0,233,29,253,253,224,239,255,0,184,168,1,0,0,0,247,247,247,247,255,0,136,168,1,0,59,38,245,230,245,208,255,0,120,168,1,0,63,129,215,161,215,106,255,0,104,168,1,0,68,197,146,77,146,33,255,0,64,168,1,0,230,220,197,197,27,125,255,0,16,168,1,0,232,118,222,222,119,174,255,0,208,167,1,0,229,62,241,241,182,218,255,0,232,166,1,0,233,29,253,253,224,239,255,0,184,166,1,0,59,38,245,230,245,208,255,0,168,166,1,0,61,103,225,184,225,134,255,0,152,166,1,0,63,166,188,127,188,65,255,0,120,166,1,0,68,197,146,77,146,33,255,0,104,166,1,0,230,220,197,197,27,125,255,0,88,166,1,0,232,118,222,222,119,174,255,0,40,166,1,0,229,62,241,241,182,218,255,0,8,166,1,0,233,29,253,253,224,239,255,0,216,165,1,0,0,0,247,247,247,247,255,0,48,165,1,0,59,38,245,230,245,208,255,0,248,164,1,0,61,103,225,184,225,134,255,0,232,164,1,0,63,166,188,127,188,65,255,0,216,164,1,0,68,197,146,77,146,33,255,0,184,164,1,0,206,255,75,64,0,75,255,0,168,164,1,0,101,255,68,0,68,27,255,0,152,164,1,0,206,173,131,118,42,131,255,0,80,164,1,0,199,87,171,153,112,171,255,0,8,164,1,0,199,51,207,194,165,207,255,0,240,163,1,0,210,21,232,231,212,232,255,0,96,163,1,0,76,30,240,217,240,211,255,0,80,163,1,0,80,68,219,166,219,160,255,0,64,163,1,0,88,123,174,90,174,97,255,0,48,163,1,0,97,197,120,27,120,55,255,0,16,163,1,0,206,255,75,64,0,75,255,0,0,163,1,0,97,197,120,27,120,55,255,0,176,162,1,0,101,255,68,0,68,27,255,0,144,162,1,0,206,173,131,118,42,131,255,0,112,162,1,0,199,87,171,153,112,171,255,0,56,162,1,0,199,51,207,194,165,207,255,0,136,161,1,0,210,21,232,231,212,232,255,0,104,161,1,0,0,0,247,247,247,247,255,0,88,161,1,0,76,30,240,217,240,211,255,0,72,161,1,0,80,68,219,166,219,160,255,0,24,161,1,0,88,123,174,90,174,97,255,0,8,161,1,0,196,70,195,175,141,195,255,0,248,160,1,0,0,0,247,247,247,247,255,0,216,160,1,0,82,90,191,127,191,123,255,0,128,160,1,0,201,168,148,123,50,148,255,0,104,160,1,0,199,51,207,194,165,207,255,0,192,159,1,0,80,68,219,166,219,160,255,0,160,159,1,0,102,255,136,0,136,55,255,0,144,159,1,0,201,168,148,123,50,148,255,0,128,159,1,0,199,51,207,194,165,207,255,0,96,159,1,0,0,0,247,247,247,247,255,0,80,159,1,0,80,68,219,166,219,160,255,0,64,159,1,0,102,255,136,0,136,55,255,0,32,159,1,0,206,173,131,118,42,131,255,0,232,158,1,0,196,70,195,175,141,195,255,0,208,158,1,0,210,21,232,231,212,232,255,0,24,158,1,0,76,30,240,217,240,211,255,0,248,157,1,0,82,90,191,127,191,123,255,0,232,157,1,0,97,197,120,27,120,55,255,0,216,157,1,0,206,173,131,118,42,131,255,0,184,157,1,0,196,70,195,175,141,195,255,0,168,157,1,0,210,21,232,231,212,232,255,0,152,157,1,0,0,0,247,247,247,247,255,0,128,157,1,0,76,30,240,217,240,211,255,0,72,157,1,0,82,90,191,127,191,123,255,0,248,156,1,0,97,197,120,27,120,55,255,0,72,156,1,0,206,173,131,118,42,131,255,0,40,156,1,0,199,87,171,153,112,171,255,0,24,156,1,0,199,51,207,194,165,207,255,0,0,156,1,0,210,21,232,231,212,232,255,0,224,155,1,0,76,30,240,217,240,211,255,0,168,155,1,0,80,68,219,166,219,160,255,0,80,155,1,0,88,123,174,90,174,97,255,0,48,155,1,0,97,197,120,27,120,55,255,0,0,155,1,0,206,173,131,118,42,131,255,0,208,154,1,0,199,87,171,153,112,171,255,0,48,154,1,0,199,51,207,194,165,207,255,0,248,153,1,0,210,21,232,231,212,232,255,0,232,153,1,0,0,0,247,247,247,247,255,0,168,153,1,0,76,30,240,217,240,211,255,0,24,153,1,0,80,68,219,166,219,160,255,0,8,153,1,0,88,123,174,90,174,97,255,0,232,152,1,0,97,197,120,27,120,55,255,0,208,152,1,0,189,11,242,236,231,242,255,0,184,152,1,0,151,61,219,166,189,219,255,0,112,152,1,0,141,197,190,43,140,190,255,0,152,151,1,0,185,8,246,241,238,246,255,0,88,151,1,0,155,40,225,189,201,225,255,0,40,151,1,0,145,112,207,116,169,207,255,0,8,151,1,0,143,247,176,5,112,176,255,0,200,150,1,0,185,8,246,241,238,246,255,0,184,150,1,0,155,40,225,189,201,225,255,0,168,150,1,0,145,112,207,116,169,207,255,0,112,150,1,0,141,197,190,43,140,190,255,0,80,150,1,0,143,247,141,4,90,141,255,0,232,149,1,0,185,8,246,241,238,246,255,0,24,149,1,0,168,24,230,208,209,230,255,0,248,148,1,0,151,61,219,166,189,219,255,0,168,148,1,0,145,112,207,116,169,207,255,0,152,148,1,0,141,197,190,43,140,190,255,0,112,148,1,0,143,247,141,4,90,141,255,0,96,148,1,0,185,8,246,241,238,246,255,0,72,148,1,0,168,24,230,208,209,230,255,0,48,148,1,0,151,61,219,166,189,219,255,0,248,147,1,0,145,112,207,116,169,207,255,0,200,147,1,0,142,183,192,54,144,192,255,0,64,147,1,0,143,247,176,5,112,176,255,0,32,147,1,0,143,248,123,3,78,123,255,0,16,147,1,0,233,8,255,255,247,251,255,0,0,147,1,0,189,11,242,236,231,242,255,0,224,146,1,0,168,24,230,208,209,230,255,0,208,146,1,0,151,61,219,166,189,219,255,0,192,146,1,0,145,112,207,116,169,207,255,0,160,146,1,0,142,183,192,54,144,192,255,0,136,146,1,0,143,247,176,5,112,176,255,0,40,146,1,0,143,248,123,3,78,123,255,0,104,145,1,0,233,8,255,255,247,251,255,0,72,145,1,0,189,11,242,236,231,242,255,0,56,145,1,0,168,24,230,208,209,230,255,0,40,145,1,0,151,61,219,166,189,219,255,0,8,145,1,0,145,112,207,116,169,207,255,0,248,144,1,0,142,183,192,54,144,192,255,0,232,144,1,0,143,247,176,5,112,176,255,0,192,144,1,0,143,247,141,4,90,141,255,0,160,144,1,0,143,249,88,2,56,88,255,0,144,144,1,0,200,14,240,236,226,240,255,0,16,144,1,0,151,61,219,166,189,219,255,0,240,143,1,0,130,208,153,28,144,153,255,0,224,143,1,0,207,8,247,246,239,247,255,0,208,143,1,0,155,40,225,189,201,225,255,0,176,143,1,0,143,128,207,103,169,207,255,0,160,143,1,0,130,251,138,2,129,138,255,0,144,143,1,0,207,8,247,246,239,247,255,0,88,143,1,0,155,40,225,189,201,225,255,0,72,143,1,0,143,128,207,103,169,207,255,0,56,143,1,0,130,208,153,28,144,153,255,0,168,142,1,0,119,252,108,1,108,89,255,0,136,142,1,0,207,8,247,246,239,247,255,0,120,142,1,0,168,24,230,208,209,230,255,0,104,142,1,0,151,61,219,166,189,219,255,0,56,142,1,0,143,128,207,103,169,207,255,0,40,142,1,0,130,208,153,28,144,153,255,0,24,142,1,0,119,252,108,1,108,89,255,0,248,141,1,0,207,8,247,246,239,247,255,0,216,141,1,0,168,24,230,208,209,230,255,0,200,141,1,0,151,61,219,166,189,219,255,0,72,141,1,0,143,128,207,103,169,207,255,0,40,141,1,0,142,183,192,54,144,192,255,0,24,141,1,0,130,251,138,2,129,138,255,0,8,141,1,0,118,252,100,1,100,80,255,0,232,140,1,0,233,8,255,255,247,251,255,0,216,140,1,0,200,14,240,236,226,240,255,0,200,140,1,0,168,24,230,208,209,230,255,0,168,140,1,0,151,61,219,166,189,219,255,0,136,140,1,0,143,128,207,103,169,207,255,0,120,140,1,0,142,183,192,54,144,192,255,0,208,139,1,0,130,251,138,2,129,138,255,0,192,139,1,0,118,252,100,1,100,80,255,0,176,139,1,0,233,8,255,255,247,251,255,0,160,139,1,0,200,14,240,236,226,240,255,0,128,139,1,0,168,24,230,208,209,230,255,0,112,139,1,0,151,61,219,166,189,219,255,0,96,139,1,0,143,128,207,103,169,207,255,0,64,139,1,0,142,183,192,54,144,192,255,0,8,139,1,0,130,251,138,2,129,138,255,0,248,138,1,0,119,252,108,1,108,89,255,0,56,138,1,0,117,251,70,1,70,54,255,0,32,138,1,0,18,238,127,127,59,8,255,0,16,138,1,0,195,255,75,45,0,75,255,0,248,137,1,0,20,246,179,179,88,6,255,0,216,137,1,0,22,232,224,224,130,20,255,0,160,137,1,0,23,155,253,253,184,99,255,0,240,136,1,0,24,72,254,254,224,182,255,0,208,136,1,0,165,20,235,216,218,235,255,0,160,136,1,0,177,47,210,178,171,210,255,0,144,136,1,0,179,84,172,128,115,172,255,0,240,135,1,0,189,181,136,84,39,136,255,0,184,135,1,0,18,238,127,127,59,8,255,0,168,135,1,0,189,181,136,84,39,136,255,0,96,135,1,0,195,255,75,45,0,75,255,0,16,135,1,0,20,246,179,179,88,6,255,0,216,134,1,0,22,232,224,224,130,20,255,0,192,134,1,0,23,155,253,253,184,99,255,0,152,134,1,0,24,72,254,254,224,182,255,0,104,134,1,0,0,0,247,247,247,247,255,0,72,134,1,0,165,20,235,216,218,235,255,0,200,133,1,0,177,47,210,178,171,210,255,0,168,133,1,0,179,84,172,128,115,172,255,0,120,133,1,0,23,187,241,241,163,64,255,0,88,133,1,0,0,0,247,247,247,247,255,0,248,132,1,0,178,69,195,153,142,195,255,0,232,132,1,0,17,253,230,230,97,1,255,0,216,132,1,0,23,155,253,253,184,99,255,0,168,132,1,0,177,47,210,178,171,210,255,0,128,132,1,0,185,155,153,94,60,153,255,0,104,132,1,0,17,253,230,230,97,1,255,0,8,132,1,0,23,155,253,253,184,99,255,0,224,131,1,0,0,0,247,247,247,247,255,0,192,131,1,0,177,47,210,178,171,210,255,0,176,131,1,0,185,155,153,94,60,153,255,0,144,131,1,0,20,246,179,179,88,6,255,0,128,131,1,0,23,187,241,241,163,64,255,0,104,131,1,0,24,72,254,254,224,182,255,0,80,131,1,0,165,20,235,216,218,235,255,0,48,131,1,0,178,69,195,153,142,195,255,0,16,131,1,0,189,181,136,84,39,136,255,0,56,130,1,0,20,246,179,179,88,6,255,0,16,130,1,0,23,187,241,241,163,64,255,0,0,130,1,0,24,72,254,254,224,182,255,0,240,129,1,0,0,0,247,247,247,247,255,0,208,129,1,0,165,20,235,216,218,235,255,0,192,129,1,0,178,69,195,153,142,195,255,0,176,129,1,0,189,181,136,84,39,136,255,0,136,129,1,0,20,246,179,179,88,6,255,0,104,129,1,0,22,232,224,224,130,20,255,0,88,129,1,0,23,155,253,253,184,99,255,0,208,128,1,0,24,72,254,254,224,182,255,0,184,128,1,0,165,20,235,216,218,235,255,0,168,128,1,0,177,47,210,178,171,210,255,0,152,128,1,0,179,84,172,128,115,172,255,0,120,128,1,0,189,181,136,84,39,136,255,0,104,128,1,0,20,246,179,179,88,6,255,0,88,128,1,0,22,232,224,224,130,20,255,0,64,128,1,0,23,155,253,253,184,99,255,0,24,128,1,0,24,72,254,254,224,182,255,0,8,128,1,0,0,0,247,247,247,247,255,0,96,127,1,0,165,20,235,216,218,235,255,0,72,127,1,0,177,47,210,178,171,210,255,0,56,127,1,0,179,84,172,128,115,172,255,0,40,127,1,0,189,181,136,84,39,136,255,0,8,127,1,0,188,14,239,231,225,239,255,0,248,126,1,0,214,67,201,201,148,199,255,0,232,126,1,0,234,222,221,221,28,119,255,0,192,126,1,0,185,8,246,241,238,246,255,0,152,126,1,0,211,41,216,215,181,216,255,0,136,126,1,0,228,139,223,223,101,176,255,0,216,125,1,0,239,232,206,206,18,86,255,0,192,125,1,0,185,8,246,241,238,246,255,0,176,125,1,0,211,41,216,215,181,216,255,0,160,125,1,0,228,139,223,223,101,176,255,0,112,125,1,0,234,222,221,221,28,119,255,0,96,125,1,0,236,255,152,152,0,67,255,0,80,125,1,0,185,8,246,241,238,246,255,0,32,125,1,0,204,38,218,212,185,218,255,0,0,125,1,0,214,67,201,201,148,199,255,0,240,124,1,0,228,139,223,223,101,176,255,0,40,124,1,0,234,222,221,221,28,119,255,0,8,124,1,0,236,255,152,152,0,67,255,0,248,123,1,0,185,8,246,241,238,246,255,0,232,123,1,0,204,38,218,212,185,218,255,0,200,123,1,0,214,67,201,201,148,199,255,0,184,123,1,0,228,139,223,223,101,176,255,0,168,123,1,0,233,209,231,231,41,138,255,0,128,123,1,0,239,232,206,206,18,86,255,0,88,123,1,0,236,255,145,145,0,63,255,0,72,123,1,0,195,5,249,247,244,249,255,0,128,122,1,0,188,14,239,231,225,239,255,0,104,122,1,0,204,38,218,212,185,218,255,0,88,122,1,0,214,67,201,201,148,199,255,0,72,122,1,0,228,139,223,223,101,176,255,0,40,122,1,0,233,209,231,231,41,138,255,0,24,122,1,0,239,232,206,206,18,86,255,0,8,122,1,0,236,255,145,145,0,63,255,0,216,121,1,0,195,5,249,247,244,249,255,0,160,121,1,0,188,14,239,231,225,239,255,0,144,121,1,0,204,38,218,212,185,218,255,0,32,121,1,0,214,67,201,201,148,199,255,0,240,120,1,0,228,139,223,223,101,176,255,0,224,120,1,0,233,209,231,231,41,138,255,0,200,120,1,0,239,232,206,206,18,86,255,0,168,120,1,0,236,255,152,152,0,67,255,0,128,120,1,0,242,255,103,103,0,31,255,0,56,120,1,0,180,8,245,239,237,245,255,0,8,120,1,0,168,37,220,188,189,220,255,0,216,119,1,0,176,100,177,117,107,177,255,0,200,119,1,0,182,7,247,242,240,247,255,0,120,119,1,0,173,28,226,203,201,226,255,0,64,119,1,0,173,58,200,158,154,200,255,0,48,119,1,0,182,128,163,106,81,163,255,0,232,118,1,0,182,7,247,242,240,247,255,0,176,118,1,0,173,28,226,203,201,226,255,0,144,118,1,0,173,58,200,158,154,200,255,0,128,118,1,0,176,100,177,117,107,177,255,0,80,118,1,0,188,185,143,84,39,143,255,0,8,118,1,0,182,7,247,242,240,247,255,0,248,117,1,0,170,18,235,218,218,235,255,0,144,117,1,0,168,37,220,188,189,220,255,0,96,117,1,0,173,58,200,158,154,200,255,0,48,117,1,0,176,100,177,117,107,177,255,0,16,117,1,0,188,185,143,84,39,143,255,0,224,116,1,0,182,7,247,242,240,247,255,0,200,116,1,0,170,18,235,218,218,235,255,0,168,116,1,0,168,37,220,188,189,220,255,0,112,116,1,0,173,58,200,158,154,200,255,0,64,116,1,0,172,83,186,128,125,186,255,0,48,116,1,0,182,128,163,106,81,163,255,0,152,115,1,0,190,216,134,74,20,134,255,0,128,115,1,0,191,2,253,252,251,253,255,0,104,115,1,0,180,8,245,239,237,245,255,0,88,115,1,0,170,18,235,218,218,235,255,0,8,115,1,0,168,37,220,188,189,220,255,0,248,114,1,0,173,58,200,158,154,200,255,0,232,114,1,0,172,83,186,128,125,186,255,0,176,114,1,0,182,128,163,106,81,163,255,0,144,114,1,0,190,216,134,74,20,134,255,0,112,114,1,0,191,2,253,252,251,253,255,0,32,114,1,0,180,8,245,239,237,245,255,0,8,114,1,0,170,18,235,218,218,235,255,0,248,113,1,0,168,37,220,188,189,220,255,0,232,113,1,0,173,58,200,158,154,200,255,0,200,113,1,0,172,83,186,128,125,186,255,0,184,113,1,0,182,128,163,106,81,163,255,0,168,113,1,0,188,185,143,84,39,143,255,0,128,113,1,0,191,255,125,63,0,125,255,0,88,113,1,0,242,255,103,103,0,31,255,0,72,113,1,0,150,241,97,5,48,97,255,0,224,112,1,0,249,220,178,178,24,43,255,0,200,112,1,0,5,163,214,214,96,77,255,0,184,112,1,0,13,119,244,244,165,130,255,0,168,112,1,0,15,54,253,253,219,199,255,0,136,112,1,0,142,32,240,209,229,240,255,0,120,112,1,0,141,87,222,146,197,222,255,0,104,112,1,0,143,167,195,67,147,195,255,0,72,112,1,0,148,206,172,33,102,172,255,0,48,112,1,0,242,255,103,103,0,31,255,0,32,112,1,0,148,206,172,33,102,172,255,0,208,111,1,0,150,241,97,5,48,97,255,0,184,111,1,0,249,220,178,178,24,43,255,0,168,111,1,0,5,163,214,214,96,77,255,0,152,111,1,0,13,119,244,244,165,130,255,0,120,111,1,0,15,54,253,253,219,199,255,0,104,111,1,0,0,0,247,247,247,247,255,0,88,111,1,0,142,32,240,209,229,240,255,0,40,111,1,0,141,87,222,146,197,222,255,0,16,111,1,0,143,167,195,67,147,195,255,0,224,110,1,0,12,150,239,239,138,98,255,0,168,110,1,0,0,0,247,247,247,247,255,0,144,110,1,0,143,128,207,103,169,207,255,0,128,110,1,0,248,255,202,202,0,32,255,0,112,110,1,0,13,119,244,244,165,130,255,0,192,173,2,0,141,87,222,146,197,222,255,0,176,173,2,0,143,247,176,5,113,176,255,0,160,173,2,0,248,255,202,202,0,32,255,0,112,173,2,0,13,119,244,244,165,130,255,0,40,173,2,0,0,0,247,247,247,247,255,0,24,173,2,0,141,87,222,146,197,222,255,0,160,172,2,0,143,247,176,5,113,176,255,0,136,172,2,0,249,220,178,178,24,43,255,0,120,172,2,0,12,150,239,239,138,98,255,0,104,172,2,0,15,54,253,253,219,199,255,0,72,172,2,0,142,32,240,209,229,240,255,0,56,172,2,0,143,128,207,103,169,207,255,0,40,172,2,0,148,206,172,33,102,172,255,0,248,171,2,0,249,220,178,178,24,43,255,0,216,171,2,0,12,150,239,239,138,98,255,0,200,171,2,0,15,54,253,253,219,199,255,0,144,171,2,0,0,0,247,247,247,247,255,0,120,171,2,0,142,32,240,209,229,240,255,0,104,171,2,0,143,128,207,103,169,207,255,0,88,171,2,0,148,206,172,33,102,172,255,0,56,171,2,0,249,220,178,178,24,43,255,0,40,171,2,0,5,163,214,214,96,77,255,0,24,171,2,0,13,119,244,244,165,130,255,0,232,170,2,0,15,54,253,253,219,199,255,0,160,170,2,0,142,32,240,209,229,240,255,0,144,170,2,0,141,87,222,146,197,222,255,0,32,170,2,0,143,167,195,67,147,195,255,0,0,170,2,0,148,206,172,33,102,172,255,0,240,169,2,0,249,220,178,178,24,43,255,0,224,169,2,0,5,163,214,214,96,77,255,0,192,169,2,0,13,119,244,244,165,130,255,0,160,169,2,0,15,54,253,253,219,199,255,0,40,169,2,0,0,0,247,247,247,247,255,0,0,169,2,0,142,32,240,209,229,240,255,0,216,168,2,0,141,87,222,146,197,222,255,0,200,168,2,0,143,167,195,67,147,195,255,0,96,168,2,0,148,206,172,33,102,172,255,0,48,168,2,0,242,255,103,103,0,31,255,0,32,168,2,0,0,0,26,26,26,26,255,0,176,167,2,0,249,220,178,178,24,43,255,0,128,167,2,0,5,163,214,214,96,77,255,0,64,167,2,0,13,119,244,244,165,130,255,0,48,167,2,0,15,54,253,253,219,199,255,0,8,167,2,0,0,0,224,224,224,224,255,0,184,166,2,0,0,0,186,186,186,186,255,0,168,166,2,0,0,0,135,135,135,135,255,0,48,166,2,0,0,0,77,77,77,77,255,0,216,165,2,0,242,255,103,103,0,31,255,0,184,165,2,0,0,0,77,77,77,77,255,0,144,165,2,0,0,0,26,26,26,26,255,0,104,165,2,0,249,220,178,178,24,43,255,0,72,165,2,0,5,163,214,214,96,77,255,0,56,165,2,0,13,119,244,244,165,130,255,0,16,165,2,0,15,54,253,253,219,199,255,0,248,164,2,0,0,0,255,255,255,255,255,0,208,164,2,0,0,0,224,224,224,224,255,0,88,164,2,0,0,0,186,186,186,186,255,0,72,164,2,0,0,0,135,135,135,135,255,0,48,164,2,0,12,150,239,239,138,98,255,0,32,164,2,0,0,0,255,255,255,255,255,0,0,164,2,0,0,0,153,153,153,153,255,0,240,163,2,0,248,255,202,202,0,32,255,0,224,163,2,0,13,119,244,244,165,130,255,0,184,163,2,0,0,0,186,186,186,186,255,0,168,163,2,0,0,0,64,64,64,64,255,0,136,163,2,0,248,255,202,202,0,32,255,0,208,162,2,0,13,119,244,244,165,130,255,0,144,162,2,0,0,0,255,255,255,255,255,0,128,162,2,0,0,0,186,186,186,186,255,0,112,162,2,0,0,0,64,64,64,64,255,0,80,162,2,0,249,220,178,178,24,43,255,0,64,162,2,0,12,150,239,239,138,98,255,0,48,162,2,0,15,54,253,253,219,199,255,0,0,162,2,0,0,0,224,224,224,224,255,0,232,161,2,0,0,0,153,153,153,153,255,0,216,161,2,0,0,0,77,77,77,77,255,0,96,161,2,0,249,220,178,178,24,43,255,0,48,161,2,0,12,150,239,239,138,98,255,0,32,161,2,0,15,54,253,253,219,199,255,0,16,161,2,0,0,0,255,255,255,255,255,0,240,160,2,0,0,0,224,224,224,224,255,0,224,160,2,0,0,0,153,153,153,153,255,0,208,160,2,0,0,0,77,77,77,77,255,0,192,160,2,0,249,220,178,178,24,43,255,0,168,160,2,0,5,163,214,214,96,77,255,0,152,160,2,0,13,119,244,244,165,130,255,0,16,160,2,0,15,54,253,253,219,199,255,0,224,159,2,0,0,0,224,224,224,224,255,0,208,159,2,0,0,0,186,186,186,186,255,0,192,159,2,0,0,0,135,135,135,135,255,0,160,159,2,0,0,0,77,77,77,77,255,0,144,159,2,0,249,220,178,178,24,43,255,0,128,159,2,0,5,163,214,214,96,77,255,0,112,159,2,0,13,119,244,244,165,130,255,0,80,159,2,0,15,54,253,253,219,199,255,0,64,159,2,0,0,0,255,255,255,255,255,0,144,158,2,0,0,0,224,224,224,224,255,0,128,158,2,0,0,0,186,186,186,186,255,0,112,158,2,0,0,0,135,135,135,135,255,0,96,158,2,0,0,0,77,77,77,77,255,0,48,158,2,0,3,32,253,253,224,221,255,0,32,158,2,0,244,92,250,250,159,181,255,0,16,158,2,0,227,220,197,197,27,138,255,0,0,158,2,0,13,28,254,254,235,226,255,0,232,157,2,0,252,72,251,251,180,185,255,0,216,157,2,0,238,147,247,247,104,161,255,0,144,157,2,0,224,253,174,174,1,126,255,0,88,157,2,0,13,28,254,254,235,226,255,0,72,157,2,0,252,72,251,251,180,185,255,0,56,157,2,0,238,147,247,247,104,161,255,0,24,157,2,0,227,220,197,197,27,138,255,0,8,157,2,0,213,252,122,122,1,119,255,0,248,156,2,0,13,28,254,254,235,226,255,0,232,156,2,0,3,60,252,252,197,192,255,0,200,156,2,0,244,92,250,250,159,181,255,0,184,156,2,0,238,147,247,247,104,161,255,0,88,156,2,0,227,220,197,197,27,138,255,0,40,156,2,0,213,252,122,122,1,119,255,0,24,156,2,0,13,28,254,254,235,226,255,0,8,156,2,0,3,60,252,252,197,192,255,0,232,155,2,0,244,92,250,250,159,181,255,0,216,155,2,0,238,147,247,247,104,161,255,0,200,155,2,0,230,195,221,221,52,151,255,0,184,155,2,0,224,253,174,174,1,126,255,0,104,155,2,0,213,252,122,122,1,119,255,0,88,155,2,0,14,12,255,255,247,243,255,0,0,155,2,0,3,32,253,253,224,221,255,0,184,154,2,0,3,60,252,252,197,192,255,0,168,154,2,0,244,92,250,250,159,181,255,0,136,154,2,0,238,147,247,247,104,161,255,0,104,154,2,0,230,195,221,221,52,151,255,0,88,154,2,0,224,253,174,174,1,126,255,0,48,154,2,0,213,252,122,122,1,119,255,0,16,154,2,0,14,12,255,255,247,243,255,0,240,153,2,0,3,32,253,253,224,221,255,0,224,153,2,0,3,60,252,252,197,192,255,0,104,153,2,0,244,92,250,250,159,181,255,0,48,153,2,0,238,147,247,247,104,161,255,0,32,153,2,0,230,195,221,221,52,151,255,0,248,152,2,0,224,253,174,174,1,126,255,0,168,152,2,0,213,252,122,122,1,119,255,0,136,152,2,0,199,255,106,73,0,106,255,0,120,152,2,0,245,255,165,165,0,38,255,0,104,152,2,0,167,171,149,49,54,149,255,0,72,152,2,0,2,208,215,215,48,39,255,0,56,152,2,0,10,184,244,244,109,67,255,0,200,151,2,0,20,157,253,253,174,97,255,0,168,151,2,0,30,110,254,254,224,144,255,0,136,151,2,0,136,24,248,224,243,248,255,0,56,151,2,0,138,67,233,171,217,233,255,0,16,151,2,0,143,113,209,116,173,209,255,0,248,150,2,0,151,157,180,69,117,180,255,0,176,150,2,0,245,255,165,165,0,38,255,0,96,150,2,0,151,157,180,69,117,180,255,0,80,150,2,0,167,171,149,49,54,149,255,0,64,150,2,0,2,208,215,215,48,39,255,0,200,149,2,0,10,184,244,244,109,67,255,0,168,149,2,0,20,157,253,253,174,97,255,0,136,149,2,0,30,110,254,254,224,144,255,0,120,149,2,0,42,64,255,255,255,191,255,0,80,149,2,0,136,24,248,224,243,248,255,0,64,149,2,0,138,67,233,171,217,233,255,0,48,149,2,0,143,113,209,116,173,209,255,0,224,148,2,0,13,164,252,252,141,89,255,0,200,148,2,0,42,64,255,255,255,191,255,0,168,148,2,0,143,86,219,145,191,219,255,0,88,148,2,0,254,225,215,215,25,28,255,0,56,148,2,0,20,157,253,253,174,97,255,0,40,148,2,0,138,67,233,171,217,233,255,0,24,148,2,0,145,193,182,44,123,182,255,0,248,147,2,0,254,225,215,215,25,28,255,0,232,147,2,0,20,157,253,253,174,97,255,0,216,147,2,0,42,64,255,255,255,191,255,0,152,147,2,0,138,67,233,171,217,233,255,0,128,147,2,0,145,193,182,44,123,182,255,0,112,147,2,0,2,208,215,215,48,39,255,0,32,147,2,0,13,164,252,252,141,89,255,0,0,147,2,0,30,110,254,254,224,144,255,0,240,146,2,0,136,24,248,224,243,248,255,0,200,146,2,0,143,86,219,145,191,219,255,0,168,146,2,0,151,157,180,69,117,180,255,0,152,146,2,0,2,208,215,215,48,39,255,0,136,146,2,0,13,164,252,252,141,89,255,0,104,146,2,0,30,110,254,254,224,144,255,0,80,146,2,0,42,64,255,255,255,191,255,0,64,146,2,0,136,24,248,224,243,248,255,0,216,145,2,0,143,86,219,145,191,219,255,0,200,145,2,0,151,157,180,69,117,180,255,0,184,145,2,0,2,208,215,215,48,39,255,0,168,145,2,0,10,184,244,244,109,67,255,0,136,145,2,0,20,157,253,253,174,97,255,0,120,145,2,0,30,110,254,254,224,144,255,0,104,145,2,0,136,24,248,224,243,248,255,0,40,145,2,0,138,67,233,171,217,233,255,0,16,145,2,0,143,113,209,116,173,209,255,0,160,144,2,0,151,157,180,69,117,180,255,0,88,144,2,0,2,208,215,215,48,39,255,0,16,144,2,0,10,184,244,244,109,67,255,0,0,144,2,0,20,157,253,253,174,97,255,0,240,143,2,0,30,110,254,254,224,144,255,0,192,143,2,0,42,64,255,255,255,191,255,0,176,143,2,0,136,24,248,224,243,248,255,0,160,143,2,0,138,67,233,171,217,233,255,0,88,143,2,0,143,113,209,116,173,209,255,0,64,143,2,0,151,157,180,69,117,180,255,0,48,143,2,0,245,255,165,165,0,38,255,0,224,142,2,0,107,255,104,0,104,55,255,0,192,142,2,0,2,208,215,215,48,39,255,0,176,142,2,0,10,184,244,244,109,67,255,0,160,142,2,0,20,157,253,253,174,97,255,0,128,142,2,0,31,115,254,254,224,139,255,0,112,142,2,0,51,106,239,217,239,139,255,0,96,142,2,0,62,130,217,166,217,106,255,0,56,142,2,0,83,121,189,102,189,99,255,0,32,142,2,0,103,211,152,26,152,80,255,0,16,142,2,0,245,255,165,165,0,38,255,0,208,141,2,0,103,211,152,26,152,80,255,0,184,141,2,0,107,255,104,0,104,55,255,0,168,141,2,0,2,208,215,215,48,39,255,0,152,141,2,0,10,184,244,244,109,67,255,0,120,141,2,0,20,157,253,253,174,97,255,0,104,141,2,0,31,115,254,254,224,139,255,0,88,141,2,0,42,64,255,255,255,191,255,0,64,141,2,0,51,106,239,217,239,139,255,0,240,140,2,0,62,130,217,166,217,106,255,0,224,140,2,0,83,121,189,102,189,99,255,0,160,140,2,0,13,164,252,252,141,89,255,0,136,140,2,0,42,64,255,255,255,191,255,0,120,140,2,0,66,136,207,145,207,96,255,0,104,140,2,0,254,225,215,215,25,28,255,0,72,140,2,0,20,157,253,253,174,97,255,0,48,140,2,0,62,130,217,166,217,106,255,0,16,140,2,0,98,210,150,26,150,65,255,0,248,139,2,0,254,225,215,215,25,28,255,0,224,139,2,0,20,157,253,253,174,97,255,0,168,139,2,0,42,64,255,255,255,191,255,0,112,139,2,0,62,130,217,166,217,106,255,0,64,139,2,0,98,210,150,26,150,65,255,0,48,139,2,0,2,208,215,215,48,39,255,0,248,138,2,0,13,164,252,252,141,89,255,0,88,138,2,0,31,115,254,254,224,139,255,0,40,138,2,0,51,106,239,217,239,139,255,0,24,138,2,0,66,136,207,145,207,96,255,0,0,138,2,0,103,211,152,26,152,80,255,0,216,137,2,0,2,208,215,215,48,39,255,0,200,137,2,0,13,164,252,252,141,89,255,0,144,137,2,0,31,115,254,254,224,139,255,0,112,137,2,0,42,64,255,255,255,191,255,0,80,137,2,0,51,106,239,217,239,139,255,0,64,137,2,0,66,136,207,145,207,96,255,0,168,136,2,0,103,211,152,26,152,80,255,0,144,136,2,0,2,208,215,215,48,39,255,0,128,136,2,0,10,184,244,244,109,67,255,0,104,136,2,0,20,157,253,253,174,97,255,0,80,136,2,0,31,115,254,254,224,139,255,0,40,136,2,0,51,106,239,217,239,139,255,0,240,135,2,0,62,130,217,166,217,106,255,0,208,135,2,0,83,121,189,102,189,99,255,0,184,135,2,0,103,211,152,26,152,80,255,0,168,135,2,0,2,208,215,215,48,39,255,0,128,135,2,0,10,184,244,244,109,67,255,0,112,135,2,0,20,157,253,253,174,97,255,0,96,135,2,0,31,115,254,254,224,139,255,0,64,135,2,0,42,64,255,255,255,191,255,0,40,135,2,0,51,106,239,217,239,139,255,0,8,135,2,0,62,130,217,166,217,106,255,0,152,134,2,0,83,121,189,102,189,99,255,0,128,134,2,0,103,211,152,26,152,80,255,0,112,134,2,0,13,44,254,254,224,210,255,0,96,134,2,0,9,139,252,252,146,114,255,0,72,134,2,0,1,211,222,222,45,38,255,0,56,134,2,0,13,37,254,254,229,217,255,0,40,134,2,0,11,108,252,252,174,145,255,0,248,133,2,0,7,179,251,251,106,74,255,0,224,133,2,0,253,224,203,203,24,29,255,0,208,133,2,0,13,37,254,254,229,217,255,0,8,133,2,0,11,108,252,252,174,145,255,0,232,132,2,0,7,179,251,251,106,74,255,0,216,132,2,0,1,211,222,222,45,38,255,0,136,132,2,0,253,231,165,165,15,21,255,0,104,132,2,0,13,37,254,254,229,217,255,0,88,132,2,0,12,92,252,252,187,161,255,0,72,132,2,0,9,139,252,252,146,114,255,0,32,132,2,0,7,179,251,251,106,74,255,0,8,132,2,0,1,211,222,222,45,38,255,0,248,131,2,0,253,231,165,165,15,21,255,0,136,131,2,0,13,37,254,254,229,217,255,0,120,131,2,0,12,92,252,252,187,161,255,0,104,131,2,0,9,139,252,252,146,114,255,0,88,131,2,0,7,179,251,251,106,74,255,0,48,131,2,0,3,208,239,239,59,44,255,0,32,131,2,0,253,224,203,203,24,29,255,0,16,131,2,0,251,255,153,153,0,13,255,0,248,130,2,0,14,15,255,255,245,240,255,0,224,130,2,0,13,44,254,254,224,210,255,0,208,130,2,0,12,92,252,252,187,161,255,0,128,130,2,0,9,139,252,252,146,114,255,0,104,130,2,0,7,179,251,251,106,74,255,0,88,130,2,0,3,208,239,239,59,44,255,0,72,130,2,0,253,224,203,203,24,29,255,0,24,130,2,0,251,255,153,153,0,13,255,0,8,130,2,0,14,15,255,255,245,240,255,0,224,129,2,0,13,44,254,254,224,210,255,0,200,129,2,0,12,92,252,252,187,161,255,0,176,129,2,0,9,139,252,252,146,114,255,0,160,129,2,0,7,179,251,251,106,74,255,0,32,129,2,0,3,208,239,239,59,44,255,0,8,129,2,0,253,224,203,203,24,29,255,0,248,128,2,0,253,231,165,165,15,21,255,0,232,128,2,0,249,255,103,103,0,13,255,0,200,128,2,0,254,225,228,228,26,28,255,0,184,128,2,0,146,178,184,55,126,184,255,0,168,128,2,0,83,147,175,77,175,74,255,0,88,128,2,0,254,225,228,228,26,28,255,0,64,128,2,0,146,178,184,55,126,184,255,0,48,128,2,0,83,147,175,77,175,74,255,0,216,127,2,0,207,132,163,152,78,163,255,0,160,127,2,0,254,225,228,228,26,28,255,0,144,127,2,0,146,178,184,55,126,184,255,0,128,127,2,0,83,147,175,77,175,74,255,0,96,127,2,0,207,132,163,152,78,163,255,0,80,127,2,0,21,255,255,255,127,0,255,0,64,127,2,0,254,225,228,228,26,28,255,0,0,127,2,0,146,178,184,55,126,184,255,0,176,126,2,0,83,147,175,77,175,74,255,0,160,126,2,0,207,132,163,152,78,163,255,0,16,126,2,0,21,255,255,255,127,0,255,0,240,125,2,0,42,204,255,255,255,51,255,0,224,125,2,0,254,225,228,228,26,28,255,0,168,125,2,0,146,178,184,55,126,184,255,0,136,125,2,0,83,147,175,77,175,74,255,0,120,125,2,0,207,132,163,152,78,163,255,0,88,125,2,0,21,255,255,255,127,0,255,0,72,125,2,0,42,204,255,255,255,51,255,0,48,125,2,0,15,193,166,166,86,40,255,0,0,125,2,0,254,225,228,228,26,28,255,0,104,124,2,0,146,178,184,55,126,184,255,0,56,124,2,0,83,147,175,77,175,74,255,0,40,124,2,0,207,132,163,152,78,163,255,0,232,123,2,0,21,255,255,255,127,0,255,0,184,123,2,0,42,204,255,255,255,51,255,0,152,123,2,0,15,193,166,166,86,40,255,0,136,123,2,0,232,121,247,247,129,191,255,0,112,123,2,0,254,225,228,228,26,28,255,0,80,123,2,0,146,178,184,55,126,184,255,0,64,123,2,0,83,147,175,77,175,74,255,0,160,122,2,0,207,132,163,152,78,163,255,0,136,122,2,0,21,255,255,255,127,0,255,0,104,122,2,0,42,204,255,255,255,51,255,0,64,122,2,0,15,193,166,166,86,40,255,0,8,122,2,0,232,121,247,247,129,191,255,0,248,121,2,0,0,0,153,153,153,153,255,0,232,121,2,0,114,120,194,102,194,165,255,0,216,121,2,0,11,155,252,252,141,98,255,0,192,121,2,0,156,77,203,141,160,203,255,0,152,121,2,0,114,120,194,102,194,165,255,0,24,121,2,0,11,155,252,252,141,98,255,0,8,121,2,0,156,77,203,141,160,203,255,0,224,120,2,0,228,102,231,231,138,195,255,0,208,120,2,0,114,120,194,102,194,165,255,0,160,120,2,0,11,155,252,252,141,98,255,0,144,120,2,0,156,77,203,141,160,203,255,0,128,120,2,0,228,102,231,231,138,195,255,0,112,120,2,0,58,155,216,166,216,84,255,0])
.concat([88,120,2,0,114,120,194,102,194,165,255,0,56,120,2,0,11,155,252,252,141,98,255,0,240,119,2,0,156,77,203,141,160,203,255,0,224,119,2,0,228,102,231,231,138,195,255,0,168,119,2,0,58,155,216,166,216,84,255,0,112,119,2,0,34,208,255,255,217,47,255,0,72,119,2,0,114,120,194,102,194,165,255,0,56,119,2,0,11,155,252,252,141,98,255,0,40,119,2,0,156,77,203,141,160,203,255,0,24,119,2,0,228,102,231,231,138,195,255,0,240,118,2,0,58,155,216,166,216,84,255,0,224,118,2,0,34,208,255,255,217,47,255,0,136,118,2,0,25,90,229,229,196,148,255,0,72,118,2,0,114,120,194,102,194,165,255,0,56,118,2,0,11,155,252,252,141,98,255,0,40,118,2,0,156,77,203,141,160,203,255,0,232,117,2,0,228,102,231,231,138,195,255,0,208,117,2,0,58,155,216,166,216,84,255,0,192,117,2,0,34,208,255,255,217,47,255,0,16,117,2,0,25,90,229,229,196,148,255,0,248,116,2,0,0,0,179,179,179,179,255,0,192,116,2,0,120,84,211,141,211,199,255,0,128,116,2,0,211,82,189,188,128,189,255,0,112,116,2,0,42,76,255,255,255,179,255,0,80,116,2,0,175,37,218,190,186,218,255,0,64,116,2,0,4,139,251,251,128,114,255,0,32,116,2,0,144,100,211,128,177,211,255,0,192,115,2,0,22,156,253,253,180,98,255,0,176,115,2,0,58,134,222,179,222,105,255,0,144,115,2,0,233,47,252,252,205,229,255,0,128,115,2,0,0,0,217,217,217,217,255,0,80,115,2,0,120,84,211,141,211,199,255,0,32,115,2,0,211,82,189,188,128,189,255,0,8,115,2,0,77,41,235,204,235,197,255,0,248,114,2,0,42,76,255,255,255,179,255,0,232,114,2,0,175,37,218,190,186,218,255,0,184,114,2,0,4,139,251,251,128,114,255,0,152,114,2,0,144,100,211,128,177,211,255,0,136,114,2,0,22,156,253,253,180,98,255,0,112,114,2,0,58,134,222,179,222,105,255,0,88,114,2,0,233,47,252,252,205,229,255,0,72,114,2,0,0,0,217,217,217,217,255,0,216,113,2,0,120,84,211,141,211,199,255,0,200,113,2,0,211,82,189,188,128,189,255,0,176,113,2,0,77,41,235,204,235,197,255,0,160,113,2,0,37,144,255,255,237,111,255,0,128,113,2,0,42,76,255,255,255,179,255,0,112,113,2,0,175,37,218,190,186,218,255,0,96,113,2,0,4,139,251,251,128,114,255,0,80,113,2,0,144,100,211,128,177,211,255,0,56,113,2,0,22,156,253,253,180,98,255,0,0,113,2,0,58,134,222,179,222,105,255,0,152,112,2,0,233,47,252,252,205,229,255,0,136,112,2,0,0,0,217,217,217,217,255,0,120,112,2,0,120,84,211,141,211,199,255,0,104,112,2,0,42,76,255,255,255,179,255,0,72,112,2,0,175,37,218,190,186,218,255,0,56,112,2,0,120,84,211,141,211,199,255,0,40,112,2,0,42,76,255,255,255,179,255,0,24,112,2,0,175,37,218,190,186,218,255,0,208,111,2,0,4,139,251,251,128,114,255,0,192,111,2,0,120,84,211,141,211,199,255,0,88,111,2,0,42,76,255,255,255,179,255,0,64,111,2,0,175,37,218,190,186,218,255,0,48,111,2,0,4,139,251,251,128,114,255,0,24,111,2,0,144,100,211,128,177,211,255,0,248,110,2,0,120,84,211,141,211,199,255,0,104,110,2,0,42,76,255,255,255,179,255,0,80,110,2,0,175,37,218,190,186,218,255,0,64,110,2,0,4,139,251,251,128,114,255,0,48,110,2,0,144,100,211,128,177,211,255,0,248,109,2,0,22,156,253,253,180,98,255,0,128,109,2,0,120,84,211,141,211,199,255,0,88,109,2,0,42,76,255,255,255,179,255,0,72,109,2,0,175,37,218,190,186,218,255,0,40,109,2,0,4,139,251,251,128,114,255,0,224,108,2,0,144,100,211,128,177,211,255,0,208,108,2,0,22,156,253,253,180,98,255,0,184,108,2,0,58,134,222,179,222,105,255,0,168,108,2,0,120,84,211,141,211,199,255,0,136,108,2,0,42,76,255,255,255,179,255,0,120,108,2,0,175,37,218,190,186,218,255,0,80,108,2,0,4,139,251,251,128,114,255,0,48,108,2,0,144,100,211,128,177,211,255,0,16,108,2,0,22,156,253,253,180,98,255,0,184,107,2,0,58,134,222,179,222,105,255,0,96,107,2,0,233,47,252,252,205,229,255,0,80,107,2,0,120,84,211,141,211,199,255,0,64,107,2,0,42,76,255,255,255,179,255,0,48,107,2,0,175,37,218,190,186,218,255,0,24,107,2,0,4,139,251,251,128,114,255,0,248,106,2,0,144,100,211,128,177,211,255,0,192,106,2,0,22,156,253,253,180,98,255,0,160,106,2,0,58,134,222,179,222,105,255,0,96,106,2,0,233,47,252,252,205,229,255,0,80,106,2,0,0,0,217,217,217,217,255,0,32,106,2,0,237,253,158,158,1,66,255,0,232,105,2,0,177,130,162,94,79,162,255,0,216,105,2,0,250,180,213,213,62,79,255,0,200,105,2,0,10,184,244,244,109,67,255,0,176,105,2,0,20,157,253,253,174,97,255,0,144,105,2,0,31,115,254,254,224,139,255,0,80,105,2,0,49,96,245,230,245,152,255,0,64,105,2,0,79,65,221,171,221,164,255,0,48,105,2,0,114,120,194,102,194,165,255,0,248,104,2,0,143,187,189,50,136,189,255,0,176,104,2,0,237,253,158,158,1,66,255,0,144,104,2,0,143,187,189,50,136,189,255,0,128,104,2,0,177,130,162,94,79,162,255,0,112,104,2,0,250,180,213,213,62,79,255,0,88,104,2,0,10,184,244,244,109,67,255,0,64,104,2,0,20,157,253,253,174,97,255,0,0,104,2,0,31,115,254,254,224,139,255,0,240,103,2,0,42,64,255,255,255,191,255,0,224,103,2,0,49,96,245,230,245,152,255,0,208,103,2,0,79,65,221,171,221,164,255,0,176,103,2,0,114,120,194,102,194,165,255,0,160,103,2,0,13,164,252,252,141,89,255,0,144,103,2,0,42,64,255,255,255,191,255,0,128,103,2,0,81,77,213,153,213,148,255,0,104,103,2,0,254,225,215,215,25,28,255,0,80,103,2,0,20,157,253,253,174,97,255,0,224,102,2,0,79,65,221,171,221,164,255,0,208,102,2,0,143,196,186,43,131,186,255,0,192,102,2,0,254,225,215,215,25,28,255,0,168,102,2,0,20,157,253,253,174,97,255,0,136,102,2,0,42,64,255,255,255,191,255,0,120,102,2,0,79,65,221,171,221,164,255,0,104,102,2,0,143,196,186,43,131,186,255,0,72,102,2,0,250,180,213,213,62,79,255,0,48,102,2,0,13,164,252,252,141,89,255,0,32,102,2,0,31,115,254,254,224,139,255,0,224,101,2,0,49,96,245,230,245,152,255,0,200,101,2,0,81,77,213,153,213,148,255,0,184,101,2,0,143,187,189,50,136,189,255,0,168,101,2,0,250,180,213,213,62,79,255,0,128,101,2,0,13,164,252,252,141,89,255,0,112,101,2,0,31,115,254,254,224,139,255,0,96,101,2,0,42,64,255,255,255,191,255,0,80,101,2,0,49,96,245,230,245,152,255,0,56,101,2,0,81,77,213,153,213,148,255,0,40,101,2,0,143,187,189,50,136,189,255,0,216,100,2,0,250,180,213,213,62,79,255,0,200,100,2,0,10,184,244,244,109,67,255,0,184,100,2,0,20,157,253,253,174,97,255,0,168,100,2,0,31,115,254,254,224,139,255,0,136,100,2,0,49,96,245,230,245,152,255,0,120,100,2,0,79,65,221,171,221,164,255,0,104,100,2,0,114,120,194,102,194,165,255,0,88,100,2,0,143,187,189,50,136,189,255,0,64,100,2,0,250,180,213,213,62,79,255,0,48,100,2,0,10,184,244,244,109,67,255,0,240,99,2,0,20,157,253,253,174,97,255,0,208,99,2,0,31,115,254,254,224,139,255,0,192,99,2,0,42,64,255,255,255,191,255,0,176,99,2,0,49,96,245,230,245,152,255,0,128,99,2,0,79,65,221,171,221,164,255,0,112,99,2,0,114,120,194,102,194,165,255,0,96,99,2,0,143,187,189,50,136,189,255,0,80,99,2,0,147,15,255,240,248,255,255,0,248,98,2,0,24,35,250,250,235,215,255,0,224,98,2,0,127,255,255,0,255,255,255,0,136,98,2,0,113,128,255,127,255,212,255,0,104,98,2,0,127,15,255,240,255,255,255,0,88,98,2,0,42,26,245,245,245,220,255,0,72,98,2,0,23,58,255,255,228,196,255,0,32,98,2,0,0,0,0,0,0,0,255,0,192,97,2,0,25,49,255,255,235,205,255,0,168,97,2,0,170,255,255,0,0,255,255,0,152,97,2,0,192,206,226,138,43,226,255,0,128,97,2,0,0,190,165,165,42,42,255,0,112,97,2,0,23,99,222,222,184,135,255,0,56,97,2,0,128,103,160,95,158,160,255,0,16,97,2,0,63,255,255,127,255,0,255,0,0,97,2,0,17,218,210,210,105,30,255,0,240,96,2,0,11,175,255,255,127,80,255,0,184,96,2,0,154,147,237,100,149,237,255,0,152,96,2,0,33,34,255,255,248,220,255,0,136,96,2,0,246,231,220,220,20,60,255,0,120,96,2,0,127,255,255,0,255,255,255,0,96,96,2,0,170,255,139,0,0,139,255,0,80,96,2,0,127,255,139,0,139,139,255,0,0,96,2,0,30,239,184,184,134,11,255,0,224,95,2,0,0,0,169,169,169,169,255,0,208,95,2,0,85,255,100,0,100,0,255,0,176,95,2,0,0,0,169,169,169,169,255,0,88,95,2,0,39,110,189,189,183,107,255,0,56,95,2,0,212,255,139,139,0,139,255,0,32,95,2,0,58,142,107,85,107,47,255,0,0,95,2,0,23,255,255,255,140,0,255,0,232,94,2,0,198,192,204,153,50,204,255,0,192,94,2,0,0,255,139,139,0,0,255,0,136,94,2,0,10,121,233,233,150,122,255,0,112,94,2,0,85,61,188,143,188,143,255,0,80,94,2,0,175,143,139,72,61,139,255,0,56,94,2,0,127,103,79,47,79,79,255,0,8,94,2,0,127,103,79,47,79,79,255,0,240,93,2,0,128,255,209,0,206,209,255,0,224,93,2,0,199,255,211,148,0,211,255,0,208,93,2,0,232,235,255,255,20,147,255,0,176,93,2,0,138,255,255,0,191,255,255,0,144,93,2,0,0,0,105,105,105,105,255,0,40,93,2,0,0,0,105,105,105,105,255,0,8,93,2,0,148,225,255,30,144,255,255,0,248,92,2,0,0,206,178,178,34,34,255,0,224,92,2,0,28,15,255,255,250,240,255,0,184,92,2,0,85,192,139,34,139,34,255,0,168,92,2,0,212,255,255,255,0,255,255,0,144,92,2,0,0,0,220,220,220,220,255,0,128,92,2,0,170,7,255,248,248,255,255,0,104,92,2,0,35,255,255,255,215,0,255,0,88,92,2,0,30,217,218,218,165,32,255,0,232,91,2,0,0,0,128,128,128,128,255,0,208,91,2,0,85,255,128,0,128,0,255,0,184,91,2,0,59,208,255,173,255,47,255,0,168,91,2,0,0,0,128,128,128,128,255,0,136,91,2,0,85,15,255,240,255,240,255,0,120,91,2,0,233,150,255,255,105,180,255,0,104,91,2,0,0,140,205,205,92,92,255,0,88,91,2,0,194,255,130,75,0,130,255,0,64,91,2,0,42,15,255,255,255,240,255,0,48,91,2,0,38,106,240,240,230,140,255,0,192,90,2,0,170,20,250,230,230,250,255,0,168,90,2,0,240,15,255,255,240,245,255,0,152,90,2,0,64,255,252,124,252,0,255,0,128,90,2,0,38,49,255,255,250,205,255,0,96,90,2,0,137,63,230,173,216,230,255,0,72,90,2,0,0,119,240,240,128,128,255,0,56,90,2,0,127,31,255,224,255,255,255,0,24,90,2,0,42,40,250,250,250,210,255,0,0,90,2,0,0,0,211,211,211,211,255,0,240,89,2,0,85,100,238,144,238,144,255,0,192,89,2,0,0,0,211,211,211,211,255,0,176,89,2,0,248,73,255,255,182,193,255,0,152,89,2,0,12,132,255,255,160,122,255,0,128,89,2,0,125,209,178,32,178,170,255,0,80,89,2,0,143,117,250,135,206,250,255,0,56,89,2,0,148,56,153,119,136,153,255,0,32,89,2,0,148,56,153,119,136,153,255,0,8,89,2,0,151,52,222,176,196,222,255,0,232,88,2,0,42,31,255,255,255,224,255,0,216,88,2,0,85,255,255,0,255,0,255,0,80,88,2,0,85,192,205,50,205,50,255,0,64,88,2,0,21,20,250,250,240,230,255,0,48,88,2,0,212,255,255,255,0,255,255,0,32,88,2,0,0,255,128,128,0,0,255,0,248,87,2,0,113,128,205,102,205,170,255,0,232,87,2,0,170,255,205,0,0,205,255,0,200,87,2,0,204,152,211,186,85,211,255,0,176,87,2,0,183,124,219,147,112,219,255,0,144,87,2,0,103,169,179,60,179,113,255,0,120,87,2,0,176,143,238,123,104,238,255,0,40,87,2,0,111,255,250,0,250,154,255,0,16,87,2,0,125,167,209,72,209,204,255,0,248,86,2,0,228,228,199,199,21,133,255,0,224,86,2,0,170,198,112,25,25,112,255,0,192,86,2,0,106,9,255,245,255,250,255,0,176,86,2,0,4,30,255,255,228,225,255,0,160,86,2,0,26,73,255,255,228,181,255,0,136,86,2,0,25,81,255,255,222,173,255,0,40,86,2,0,170,255,128,0,0,128,255,0,24,86,2,0,27,23,253,253,245,230,255,0,192,85,2,0,42,255,128,128,128,0,255,0,152,85,2,0,56,192,142,107,142,35,255,0,136,85,2,0,27,255,255,255,165,0,255,0,120,85,2,0,11,255,255,255,69,0,255,0,88,85,2,0,214,123,218,218,112,214,255,0,64,85,2,0,38,72,238,238,232,170,255,0,208,84,2,0,85,100,251,152,251,152,255,0,160,84,2,0,127,67,238,175,238,238,255,0,136,84,2,0,241,124,219,219,112,147,255,0,88,84,2,0,26,41,255,255,239,213,255,0,0,84,2,0,20,70,255,255,218,185,255,0,240,83,2,0,20,176,205,205,133,63,255,0,224,83,2,0,247,63,255,255,192,203,255,0,208,83,2,0,212,70,221,221,160,221,255,0,168,83,2,0,132,59,230,176,224,230,255,0,152,83,2,0,212,255,128,128,0,128,255,0,136,83,2,0,0,255,255,255,0,0,255,0,120,83,2,0,0,61,188,188,143,143,255,0,48,83,2,0,159,181,225,65,105,225,255,0,24,83,2,0,17,220,139,139,69,19,255,0,232,82,2,0,4,138,250,250,128,114,255,0,208,82,2,0,19,154,244,244,164,96,255,0,192,82,2,0,103,170,139,46,139,87,255,0,120,82,2,0,17,16,255,255,245,238,255,0,64,82,2,0,13,183,160,160,82,45,255,0,32,82,2,0,0,0,192,192,192,192,255,0,16,82,2,0,139,108,235,135,206,235,255,0,216,81,2,0,175,143,205,106,90,205,255,0,192,81,2,0,148,56,144,112,128,144,255,0,152,81,2,0,148,56,144,112,128,144,255,0,56,81,2,0,0,5,255,255,250,250,255,0,32,81,2,0,106,255,255,0,255,127,255,0,248,80,2,0,146,155,180,70,130,180,255,0,232,80,2,0,24,84,210,210,180,140,255,0,176,80,2,0,127,255,128,0,128,128,255,0,160,80,2,0,212,29,216,216,191,216,255,0,144,80,2,0,6,184,255,255,99,71,255,0,120,80,2,0,123,182,224,64,224,208,255,0,96,80,2,0,212,115,238,238,130,238,255,0,64,80,2,0,27,68,245,245,222,179,255,0,0,80,2,0,0,0,255,255,255,255,255,0,240,79,2,0,0,0,245,245,245,245,255,0,224,79,2,0,42,255,255,255,255,0,255,0,200,79,2,0,56,192,205,154,205,50,255,0,168,79,2,0,45,67,252,247,252,185,255,0,152,79,2,0,68,91,221,173,221,142,255,0,136,79,2,0,98,178,163,49,163,84,255,0,120,79,2,0,42,50,255,255,255,204,255,0,96,79,2,0,62,85,230,194,230,153,255,0,80,79,2,0,85,100,198,120,198,121,255,0,24,79,2,0,99,187,132,35,132,67,255,0,8,79,2,0,42,50,255,255,255,204,255,0,248,78,2,0,62,85,230,194,230,153,255,0,232,78,2,0,85,100,198,120,198,121,255,0,200,78,2,0,98,178,163,49,163,84,255,0,184,78,2,0,107,255,104,0,104,55,255,0,168,78,2,0,42,50,255,255,255,204,255,0,152,78,2,0,55,81,240,217,240,163,255,0,128,78,2,0,68,91,221,173,221,142,255,0,112,78,2,0,85,100,198,120,198,121,255,0,16,78,2,0,98,178,163,49,163,84,255,0,0,78,2,0,107,255,104,0,104,55,255,0,240,77,2,0,42,50,255,255,255,204,255,0,224,77,2,0,55,81,240,217,240,163,255,0,192,77,2,0,68,91,221,173,221,142,255,0,176,77,2,0,85,100,198,120,198,121,255,0,160,77,2,0,96,158,171,65,171,93,255,0,144,77,2,0,99,187,132,35,132,67,255,0,128,77,2,0,108,255,90,0,90,50,255,0,112,77,2,0,42,25,255,255,255,229,255,0,56,77,2,0,45,67,252,247,252,185,255,0,40,77,2,0,55,81,240,217,240,163,255,0,24,77,2,0,68,91,221,173,221,142,255,0,8,77,2,0,85,100,198,120,198,121,255,0,224,76,2,0,96,158,171,65,171,93,255,0,208,76,2,0,99,187,132,35,132,67,255,0,192,76,2,0,108,255,90,0,90,50,255,0,176,76,2,0,42,25,255,255,255,229,255,0,152,76,2,0,45,67,252,247,252,185,255,0,136,76,2,0,55,81,240,217,240,163,255,0,72,76,2,0,68,91,221,173,221,142,255,0,56,76,2,0,85,100,198,120,198,121,255,0,40,76,2,0,96,158,171,65,171,93,255,0,24,76,2,0,99,187,132,35,132,67,255,0,248,75,2,0,107,255,104,0,104,55,255,0,232,75,2,0,110,255,69,0,69,41,255,0,216,75,2,0,49,73,248,237,248,177,255,0,200,75,2,0,117,97,205,127,205,187,255,0,176,75,2,0,144,194,184,44,127,184,255,0,160,75,2,0,42,50,255,255,255,204,255,0,96,75,2,0,99,66,218,161,218,180,255,0,80,75,2,0,132,170,196,65,182,196,255,0,64,75,2,0,150,203,168,34,94,168,255,0,48,75,2,0,42,50,255,255,255,204,255,0,16,75,2,0,99,66,218,161,218,180,255,0,0,75,2,0,132,170,196,65,182,196,255,0,240,74,2,0,144,194,184,44,127,184,255,0,224,74,2,0,164,191,148,37,52,148,255,0,200,74,2,0,42,50,255,255,255,204,255,0,152,74,2,0,69,58,233,199,233,180,255,0,96,74,2,0,117,97,205,127,205,187,255,0,80,74,2,0,132,170,196,65,182,196,255,0,64,74,2,0,144,194,184,44,127,184,255,0,48,74,2,0,164,191,148,37,52,148,255,0,16,74,2,0,42,50,255,255,255,204,255,0,248,73,2,0,69,58,233,199,233,180,255,0,200,73,2,0,117,97,205,127,205,187,255,0,184,73,2,0,132,170,196,65,182,196,255,0,160,73,2,0,139,216,192,29,145,192,255,0,144,73,2,0,150,203,168,34,94,168,255,0,232,72,2,0,158,231,132,12,44,132,255,0,216,72,2,0,42,38,255,255,255,217,255,0,200,72,2,0,49,73,248,237,248,177,255,0,184,72,2,0,69,58,233,199,233,180,255,0,96,72,2,0,117,97,205,127,205,187,255,0,64,72,2,0,132,170,196,65,182,196,255,0,48,72,2,0,139,216,192,29,145,192,255,0,32,72,2,0,150,203,168,34,94,168,255,0,232,71,2,0,158,231,132,12,44,132,255,0,216,71,2,0,42,38,255,255,255,217,255,0,192,70,2,0,49,73,248,237,248,177,255,0,160,70,2,0,69,58,233,199,233,180,255,0,128,70,2,0,117,97,205,127,205,187,255,0,104,70,2,0,132,170,196,65,182,196,255,0,32,70,2,0,139,216,192,29,145,192,255,0,8,70,2,0,150,203,168,34,94,168,255,0,248,69,2,0,164,191,148,37,52,148,255,0,184,69,2,0,158,231,88,8,29,88,255,0,160,69,2,0,37,66,255,255,247,188,255,0,128,69,2,0,28,175,254,254,196,79,255,0,240,65,2,0,16,238,217,217,95,14,255,0,224,65,2,0,42,42,255,255,255,212,255,0,192,65,2,0,28,112,254,254,217,142,255,0,176,65,2,0,22,213,254,254,153,41,255,0,136,65,2,0,15,252,204,204,76,2,255,0,120,65,2,0,42,42,255,255,255,212,255,0,104,65,2,0,28,112,254,254,217,142,255,0,64,65,2,0,22,213,254,254,153,41,255,0,40,65,2,0,16,238,217,217,95,14,255,0,8,65,2,0,13,248,153,153,52,4,255,0,176,64,2,0,42,42,255,255,255,212,255,0,160,64,2,0,31,109,254,254,227,145,255,0,144,64,2,0,28,175,254,254,196,79,255,0,128,64,2,0,22,213,254,254,153,41,255,0,96,64,2,0,16,238,217,217,95,14,255,0,80,64,2,0,13,248,153,153,52,4,255,0,64,64,2,0,42,42,255,255,255,212,255,0,48,64,2,0,31,109,254,254,227,145,255,0,24,64,2,0,28,175,254,254,196,79,255,0,8,64,2,0,22,213,254,254,153,41,255,0,192,63,2,0,18,233,236,236,112,20,255,0,176,63,2,0,15,252,204,204,76,2,255,0,160,63,2,0,12,247,140,140,45,4,255,0,144,63,2,0,42,25,255,255,255,229,255,0,104,63,2,0,37,66,255,255,247,188,255,0,88,63,2,0,31,109,254,254,227,145,255,0,72,63,2,0,28,175,254,254,196,79,255,0,56,63,2,0,22,213,254,254,153,41,255,0,32,63,2,0,18,233,236,236,112,20,255,0,16,63,2,0,15,252,204,204,76,2,255,0,176,62,2,0,12,247,140,140,45,4,255,0,144,62,2,0,42,25,255,255,255,229,255,0,128,62,2,0,37,66,255,255,247,188,255,0,112,62,2,0,31,109,254,254,227,145,255,0,80,62,2,0,28,175,254,254,196,79,255,0,64,62,2,0,22,213,254,254,153,41,255,0,48,62,2,0,18,233,236,236,112,20,255,0,32,62,2,0,15,252,204,204,76,2,255,0,8,62,2,0,13,248,153,153,52,4,255,0,248,61,2,0,13,240,102,102,37,6,255,0,104,61,2,0,34,95,255,255,237,160,255,0,88,61,2,0,24,178,254,254,178,76,255,0,72,61,2,0,5,221,240,240,59,32,255,0,56,61,2,0,42,77,255,255,255,178,255,0,8,61,2,0,29,162,254,254,204,92,255,0,248,60,2,0,17,194,253,253,141,60,255,0,232,60,2,0,254,225,227,227,26,28,255,0,216,60,2,0,42,77,255,255,255,178,255,0,192,60,2,0,29,162,254,254,204,92,255,0,176,60,2,0,17,194,253,253,141,60,255,0,16,60,2,0,5,221,240,240,59,32,255,0,0,60,2,0,246,255,189,189,0,38,255,0,240,59,2,0,42,77,255,255,255,178,255,0,224,59,2,0,30,136,254,254,217,118,255,0,192,59,2,0,24,178,254,254,178,76,255,0,176,59,2,0,17,194,253,253,141,60,255,0,160,59,2,0,5,221,240,240,59,32,255,0,144,59,2,0,246,255,189,189,0,38,255,0,120,59,2,0,42,77,255,255,255,178,255,0,104,59,2,0,30,136,254,254,217,118,255,0,232,58,2,0,24,178,254,254,178,76,255,0,216,58,2,0,17,194,253,253,141,60,255,0,200,58,2,0,7,212,252,252,78,42,255,0,184,58,2,0,254,225,227,227,26,28,255,0,152,58,2,0,245,255,177,177,0,38,255,0,136,58,2,0,42,50,255,255,255,204,255,0,120,58,2,0,34,95,255,255,237,160,255,0,104,58,2,0,30,136,254,254,217,118,255,0,80,58,2,0,24,178,254,254,178,76,255,0,56,58,2,0,17,194,253,253,141,60,255,0,224,57,2,0,7,212,252,252,78,42,255,0,208,57,2,0,254,225,227,227,26,28,255,0,192,57,2,0,245,255,177,177,0,38,255,0,176,57,2,0,42,50,255,255,255,204,255,0,144,57,2,0,34,95,255,255,237,160,255,0,120,57,2,0,30,136,254,254,217,118,255,0,32,57,2,0,24,178,254,254,178,76,255,0,16,57,2,0,17,194,253,253,141,60,255,0,232,56,2,0,7,212,252,252,78,42,255,0,216,56,2,0,254,225,227,227,26,28,255,0,136,56,2,0,246,255,189,189,0,38,255,0,120,56,2,0,242,255,128,128,0,38,255,0,32,3,2,0,147,15,255,240,248,255,255,0,56,234,1,0,24,35,250,250,235,215,255,0,56,56,2,0,23,36,255,255,239,219,255,0,248,55,2,0,23,36,238,238,223,204,255,0,232,55,2,0,23,36,205,205,192,176,255,0,216,55,2,0,24,34,139,139,131,120,255,0,216,189,1,0,113,128,255,127,255,212,255,0,192,55,2,0,113,128,255,127,255,212,255,0,104,55,2,0,113,128,238,118,238,198,255,0,72,55,2,0,113,128,205,102,205,170,255,0,48,55,2,0,113,128,139,69,139,116,255,0,160,169,1,0,127,15,255,240,255,255,255,0,248,54,2,0,127,15,255,240,255,255,255,0,240,54,2,0,127,15,238,224,238,238,255,0,232,54,2,0,127,14,205,193,205,205,255,0,216,54,2,0,127,14,139,131,139,139,255,0,72,150,1,0,42,26,245,245,245,220,255,0,120,132,1,0,23,58,255,255,228,196,255,0,152,54,2,0,23,58,255,255,228,196,255,0,144,54,2,0,23,58,238,238,213,183,255,0,128,54,2,0,22,58,205,205,183,158,255,0,120,54,2,0,23,58,139,139,125,107,255,0,152,149,2,0,0,0,0,0,0,0,255,0,224,164,2,0,25,49,255,255,235,205,255,0,200,135,2,0,170,255,255,0,0,255,255,0,96,54,2,0,170,255,255,0,0,255,255,0,80,54,2,0,170,255,238,0,0,238,255,0,72,54,2,0,170,255,205,0,0,205,255,0,24,54,2,0,170,255,139,0,0,139,255,0,56,136,2,0,192,206,226,138,43,226,255,0,168,121,2,0,0,190,165,165,42,42,255,0,16,54,2,0,0,191,255,255,64,64,255,0,8,54,2,0,0,191,238,238,59,59,255,0,0,54,2,0,0,191,205,205,51,51,255,0,248,53,2,0,0,190,139,139,35,35,255,0,8,107,2,0,23,99,222,222,184,135,255,0,224,53,2,0,23,100,255,255,211,155,255,0,208,53,2,0,23,99,238,238,197,145,255,0,152,53,2,0,23,99,205,205,170,125,255,0,120,53,2,0,23,99,139,139,115,85,255,0,208,94,2,0,128,103,160,95,158,160,255,0,104,53,2,0,131,103,255,152,245,255,255,0,72,53,2,0,131,102,238,142,229,238,255,0,56,53,2,0,131,103,205,122,197,205,255,0,40,53,2,0,131,102,139,83,134,139,255,0,168,81,2,0,63,255,255,127,255,0,255,0,16,53,2,0,63,255,255,127,255,0,255,0,0,53,2,0,63,255,238,118,238,0,255,0,200,52,2,0,63,255,205,102,205,0,255,0,184,52,2,0,63,255,139,69,139,0,255,0,144,69,2,0,17,218,210,210,105,30,255,0,168,52,2,0,17,219,255,255,127,36,255,0,136,52,2,0,17,219,238,238,118,33,255,0,120,52,2,0,17,218,205,205,102,29,255,0,104,52,2,0,17,220,139,139,69,19,255,0,192,54,2,0,11,175,255,255,127,80,255,0,88,52,2,0,7,169,255,255,114,86,255,0,80,52,2,0,6,169,238,238,106,80,255,0,48,52,2,0,6,169,205,205,91,69,255,0,40,52,2,0,6,168,139,139,62,47,255,0,8,47,2,0,154,147,237,100,149,237,255,0,80,39,2,0,33,34,255,255,248,220,255,0,248,51,2,0,33,34,255,255,248,220,255,0,232,51,2,0,34,35,238,238,232,205,255,0,216,51,2,0,34,34,205,205,200,177,255,0,200,51,2,0,35,34,139,139,136,120,255,0,216,32,2,0,246,231,220,220,20,60,255,0,24,25,2,0,127,255,255,0,255,255,255,0,160,51,2,0,127,255,255,0,255,255,255,0,152,51,2,0,127,255,238,0,238,238,255,0,144,51,2,0,127,255,205,0,205,205,255,0,136,51,2,0,127,255,139,0,139,139,255,0,136,2,2,0,30,239,184,184,134,11,255,0,104,51,2,0,30,240,255,255,185,15,255,0,88,51,2,0,30,240,238,238,173,14,255,0,72,51,2,0,30,240,205,205,149,12,255,0,48,51,2,0,30,240,139,139,101,8,255,0,16,254,1,0,85,255,100,0,100,0,255,0,248,249,1,0,39,110,189,189,183,107,255,0,16,245,1,0,58,142,107,85,107,47,255,0,248,50,2,0,58,143,255,202,255,112,255,0,232,50,2,0,58,143,238,188,238,104,255,0,200,50,2,0,58,143,205,162,205,90,255,0,184,50,2,0,58,143,139,110,139,61,255,0,32,243,1,0,23,255,255,255,140,0,255,0,168,50,2,0,21,255,255,255,127,0,255,0,144,50,2,0,21,255,238,238,118,0,255,0,112,50,2,0,21,255,205,205,102,0,255,0,40,50,2,0,21,255,139,139,69,0,255,0,48,240,1,0,198,192,204,153,50,204,255,0,24,50,2,0,198,193,255,191,62,255,255,0,8,50,2,0,198,192,238,178,58,238,255,0,232,49,2,0,198,192,205,154,50,205,255,0,200,49,2,0,198,192,139,104,34,139,255,0,208,233,1,0,10,121,233,233,150,122,255,0,168,231,1,0,85,61,188,143,188,143,255,0,120,49,2,0,85,62,255,193,255,193,255,0,104,49,2,0,85,62,238,180,238,180,255,0,48,49,2,0,85,62,205,155,205,155,255,0,32,49,2,0,85,62,139,105,139,105,255,0,112,229,1,0,175,143,139,72,61,139,255,0,184,226,1,0,127,103,79,47,79,79,255,0,224,48,2,0,127,104,255,151,255,255,255,0,136,48,2,0,127,103,238,141,238,238,255,0,120,48,2,0,127,104,205,121,205,205,255,0,104,48,2,0,127,104,139,82,139,139,255,0,112,224,1,0,127,103,79,47,79,79,255,0,72,222,1,0,128,255,209,0,206,209,255,0,72,220,1,0,199,255,211,148,0,211,255,0,152,218,1,0,232,235,255,255,20,147,255,0,224,47,2,0,232,235,255,255,20,147,255,0,200,47,2,0,232,235,238,238,18,137,255,0,80,47,2,0,232,235,205,205,16,118,255,0,64,47,2,0,231,236,139,139,10,80,255,0,152,216,1,0,138,255,255,0,191,255,255,0,48,47,2,0,138,255,255,0,191,255,255,0,24,47,2,0,138,255,238,0,178,238,255,0,240,46,2,0,138,255,205,0,154,205,255,0,176,46,2,0,138,255,139,0,104,139,255,0,96,212,1,0,0,0,105,105,105,105,255,0,208,209,1,0,0,0,105,105,105,105,255,0,152,207,1,0,148,225,255,30,144,255,255,0,120,46,2,0,148,225,255,30,144,255,255,0,104,46,2,0,148,225,238,28,134,238,255,0,88,46,2,0,148,225,205,24,116,205,255,0,72,46,2,0,148,225,139,16,78,139,255,0,128,205,1,0,0,206,178,178,34,34,255,0,48,46,2,0,0,207,255,255,48,48,255,0,184,45,2,0,0,207,238,238,44,44,255,0,168,45,2,0,0,207,205,205,38,38,255,0,152,45,2,0,0,207,139,139,26,26,255,0,216,203,1,0,28,15,255,255,250,240,255,0,104,202,1,0,85,192,139,34,139,34,255,0,232,198,1,0,0,0,220,220,220,220,255,0,48,197,1,0,170,7,255,248,248,255,255,0,240,193,1,0,35,255,255,255,215,0,255,0,120,45,2,0,35,255,255,255,215,0,255,0,112,45,2,0,35,255,238,238,201,0,255,0,8,45,2,0,35,255,205,205,173,0,255,0,0,45,2,0,35,255,139,139,117,0,255,0,136,191,1,0,30,217,218,218,165,32,255,0,240,44,2,0,30,218,255,255,193,37,255,0,208,44,2,0,30,218,238,238,180,34,255,0,192,44,2,0,30,218,205,205,155,29,255,0,176,44,2,0,30,218,139,139,105,20,255,0,112,106,2,0,0,0,192,192,192,192,255,0,160,44,2,0,0,0,0,0,0,0,255,0,152,44,2,0,0,0,3,3,3,3,255,0,112,220,1,0,0,0,26,26,26,26,255,0,72,44,2,0,0,0,255,255,255,255,255,0,64,44,2,0,0,0,28,28,28,28,255,0,56,44,2,0,0,0,31,31,31,31,255,0,8,44,2,0,0,0,33,33,33,33,255,0,0,44,2,0,0,0,36,36,36,36,255,0,192,218,1,0,0,0,38,38,38,38,255,0,248,43,2,0,0,0,41,41,41,41,255,0,232,43,2,0,0,0,43,43,43,43,255,0,224,43,2,0,0,0,46,46,46,46,255,0,136,43,2,0,0,0,48,48,48,48,255,0,128,43,2,0,0,0,5,5,5,5,255,0,208,216,1,0,0,0,51,51,51,51,255,0,104,43,2,0,0,0,54,54,54,54,255,0,64,43,2,0,0,0,56,56,56,56,255,0,56,43,2,0,0,0,59,59,59,59,255,0,48,43,2,0,0,0,61,61,61,61,255,0,136,212,1,0,0,0,64,64,64,64,255,0,32,43,2,0,0,0,66,66,66,66,255,0,24,43,2,0,0,0,69,69,69,69,255,0,240,42,2,0,0,0,71,71,71,71,255,0,232,42,2,0,0,0,74,74,74,74,255,0,224,42,2,0,0,0,8,8,8,8,255,0,232,209,1,0,0,0,77,77,77,77,255,0,200,42,2,0,0,0,79,79,79,79,255,0,192,42,2,0,0,0,82,82,82,82,255,0,184,42,2,0,0,0,84,84,84,84,255,0,176,42,2,0,0,0,87,87,87,87,255,0,184,207,1,0,0,0,89,89,89,89,255,0,160,42,2,0,0,0,92,92,92,92,255,0,88,42,2,0,0,0,94,94,94,94,255,0,80,42,2,0,0,0,97,97,97,97,255,0,72,42,2,0,0,0,99,99,99,99,255,0,64,42,2,0,0,0,10,10,10,10,255,0,200,205,1,0,0,0,102,102,102,102,255,0,40,42,2,0,0,0,105,105,105,105,255,0,32,42,2,0,0,0,107,107,107,107,255,0,24,42,2,0,0,0,110,110,110,110,255,0,8,42,2,0,0,0,112,112,112,112,255,0,0,204,1,0,0,0,115,115,115,115,255,0,192,41,2,0,0,0,117,117,117,117,255,0,168,41,2,0,0,0,120,120,120,120,255,0,160,41,2,0,0,0,122,122,122,122,255,0,152,41,2,0,0,0,125,125,125,125,255,0,120,41,2,0,0,0,13,13,13,13,255,0,160,202,1,0,0,0,127,127,127,127,255,0,72,41,2,0,0,0,130,130,130,130,255,0,64,41,2,0,0,0,133,133,133,133,255,0,48,41,2,0,0,0,135,135,135,135,255,0,40,41,2,0,0,0,138,138,138,138,255,0,192,200,1,0,0,0,140,140,140,140,255,0,248,40,2,0,0,0,143,143,143,143,255,0,240,40,2,0,0,0,145,145,145,145,255,0,232,40,2,0,0,0,148,148,148,148,255,0,192,40,2,0,0,0,150,150,150,150,255,0,144,40,2,0,0,0,15,15,15,15,255,0,88,199,1,0,0,0,153,153,153,153,255,0,136,40,2,0,0,0,156,156,156,156,255,0,64,40,2,0,0,0,158,158,158,158,255,0,56,40,2,0,0,0,161,161,161,161,255,0,8,40,2,0,0,0,163,163,163,163,255,0,112,197,1,0,0,0,166,166,166,166,255,0,232,39,2,0,0,0,168,168,168,168,255,0,224,39,2,0,0,0,171,171,171,171,255,0,144,39,2,0,0,0,173,173,173,173,255,0,136,39,2,0,0,0,176,176,176,176,255,0,128,39,2,0,0,0,18,18,18,18,255,0,8,194,1,0,0,0,179,179,179,179,255,0,96,39,2,0,0,0,181,181,181,181,255,0,72,39,2,0,0,0,184,184,184,184,255,0,248,38,2,0,0,0,186,186,186,186,255,0,240,38,2,0,0,0,189,189,189,189,255,0,168,191,1,0,0,0,191,191,191,191,255,0,208,38,2,0,0,0,194,194,194,194,255,0,184,38,2,0,0,0,196,196,196,196,255,0,176,38,2,0,0,0,199,199,199,199,255,0,168,38,2,0,0,0,201,201,201,201,255,0,160,38,2,0,0,0,20,20,20,20,255,0,120,189,1,0,0,0,204,204,204,204,255,0,144,38,2,0,0,0,207,207,207,207,255,0,80,38,2,0,0,0,209,209,209,209,255,0,72,38,2,0,0,0,212,212,212,212,255,0,64,38,2,0,0,0,214,214,214,214,255,0,104,187,1,0,0,0,217,217,217,217,255,0,40,38,2,0,0,0,219,219,219,219,255,0,32,38,2,0,0,0,222,222,222,222,255,0,24,38,2,0,0,0,224,224,224,224,255,0,16,38,2,0,0,0,227,227,227,227,255,0,0,38,2,0,0,0,23,23,23,23,255,0,32,186,1,0,0,0,229,229,229,229,255,0,216,37,2,0,0,0,232,232,232,232,255,0,208,37,2,0,0,0,235,235,235,235,255,0,200,37,2,0,0,0,237,237,237,237,255,0,192,37,2,0,0,0,240,240,240,240,255,0,104,184,1,0,0,0,242,242,242,242,255,0,168,37,2,0,0,0,245,245,245,245,255,0,160,37,2,0,0,0,247,247,247,247,255,0,152,37,2,0,0,0,250,250,250,250,255,0,136,37,2,0,0,0,252,252,252,252,255,0,104,94,2,0,85,255,255,0,255,0,255,0,96,37,2,0,85,255,255,0,255,0,255,0,88,37,2,0,85,255,238,0,238,0,255,0,80,37,2,0,85,255,205,0,205,0,255,0,72,37,2,0,85,255,139,0,139,0,255,0,248,185,1,0,59,208,255,173,255,47,255,0,64,184,1,0,0,0,192,192,192,192,255,0,48,37,2,0,0,0,0,0,0,0,255,0,40,37,2,0,0,0,3,3,3,3,255,0,24,37,2,0,0,0,26,26,26,26,255,0,16,37,2,0,0,0,255,255,255,255,255,0,240,36,2,0,0,0,28,28,28,28,255,0,232,36,2,0,0,0,31,31,31,31,255,0,224,36,2,0,0,0,33,33,33,33,255,0,216,36,2,0,0,0,36,36,36,36,255,0,184,36,2,0,0,0,38,38,38,38,255,0,176,36,2,0,0,0,41,41,41,41,255,0,168,36,2,0,0,0,43,43,43,43,255,0,160,36,2,0,0,0,46,46,46,46,255,0,144,36,2,0,0,0,48,48,48,48,255,0,136,36,2,0,0,0,5,5,5,5,255,0,96,36,2,0,0,0,51,51,51,51,255,0,88,36,2,0,0,0,54,54,54,54,255,0,80,36,2,0,0,0,56,56,56,56,255,0,72,36,2,0,0,0,59,59,59,59,255,0,48,36,2,0,0,0,61,61,61,61,255,0,40,36,2,0,0,0,64,64,64,64,255,0,32,36,2,0,0,0,66,66,66,66,255,0,24,36,2,0,0,0,69,69,69,69,255,0,8,36,2,0,0,0,71,71,71,71,255,0,0,36,2,0,0,0,74,74,74,74,255,0,208,35,2,0,0,0,8,8,8,8,255,0,200,35,2,0,0,0,77,77,77,77,255,0,192,35,2,0,0,0,79,79,79,79,255,0,184,35,2,0,0,0,82,82,82,82,255,0,160,35,2,0,0,0,84,84,84,84,255,0,152,35,2,0,0,0,87,87,87,87,255,0,144,35,2,0,0,0,89,89,89,89,255,0,136,35,2,0,0,0,92,92,92,92,255,0,120,35,2,0,0,0,94,94,94,94,255,0,112,35,2,0,0,0,97,97,97,97,255,0,96,35,2,0,0,0,99,99,99,99,255,0,80,35,2,0,0,0,10,10,10,10,255,0,72,35,2,0,0,0,102,102,102,102,255,0,64,35,2,0,0,0,105,105,105,105,255,0,40,35,2,0,0,0,107,107,107,107,255,0,8,35,2,0,0,0,110,110,110,110,255,0,248,34,2,0,0,0,112,112,112,112,255,0,240,34,2,0,0,0,115,115,115,115,255,0,112,34,2,0,0,0,117,117,117,117,255,0,104,34,2,0,0,0,120,120,120,120,255,0,88,34,2,0,0,0,122,122,122,122,255,0,80,34,2,0,0,0,125,125,125,125,255,0,72,34,2,0,0,0,13,13,13,13,255,0,64,34,2,0,0,0,127,127,127,127,255,0,240,33,2,0,0,0,130,130,130,130,255,0,208,33,2,0,0,0,133,133,133,133,255,0,200,33,2,0,0,0,135,135,135,135,255,0,192,33,2,0,0,0,138,138,138,138,255,0,176,33,2,0,0,0,140,140,140,140,255,0,168,33,2,0,0,0,143,143,143,143,255,0,120,33,2,0,0,0,145,145,145,145,255,0,104,33,2,0,0,0,148,148,148,148,255,0,88,33,2,0,0,0,150,150,150,150,255,0,80,33,2,0,0,0,15,15,15,15,255,0,40,33,2,0,0,0,153,153,153,153,255,0,16,33,2,0,0,0,156,156,156,156,255,0,8,33,2,0,0,0,158,158,158,158,255,0,240,32,2,0,0,0,161,161,161,161,255,0,224,32,2,0,0,0,163,163,163,163,255,0,208,32,2,0,0,0,166,166,166,166,255,0,168,32,2,0,0,0,168,168,168,168,255,0,160,32,2,0,0,0,171,171,171,171,255,0,120,32,2,0,0,0,173,173,173,173,255,0,112,32,2,0,0,0,176,176,176,176,255,0,88,32,2,0,0,0,18,18,18,18,255,0,80,32,2,0,0,0,179,179,179,179,255,0,72,32,2,0,0,0,181,181,181,181,255,0,64,32,2,0,0,0,184,184,184,184,255,0,48,32,2,0,0,0,186,186,186,186,255,0,40,32,2,0,0,0,189,189,189,189,255,0,0,32,2,0,0,0,191,191,191,191,255,0,248,31,2,0,0,0,194,194,194,194,255,0,240,31,2,0,0,0,196,196,196,196,255,0,232,31,2,0,0,0,199,199,199,199,255,0,208,31,2,0,0,0,201,201,201,201,255,0,200,31,2,0,0,0,20,20,20,20,255,0,192,31,2,0,0,0,204,204,204,204,255,0,184,31,2,0,0,0,207,207,207,207,255,0,168,31,2,0,0,0,209,209,209,209,255,0,160,31,2,0,0,0,212,212,212,212,255,0,136,31,2,0,0,0,214,214,214,214,255,0,128,31,2,0,0,0,217,217,217,217,255,0,120,31,2,0,0,0,219,219,219,219,255,0,112,31,2,0,0,0,222,222,222,222,255,0,64,31,2,0,0,0,224,224,224,224,255,0,56,31,2,0,0,0,227,227,227,227,255,0,48,31,2,0,0,0,23,23,23,23,255,0,40,31,2,0,0,0,229,229,229,229,255,0,24,31,2,0,0,0,232,232,232,232,255,0,16,31,2,0,0,0,235,235,235,235,255,0,248,30,2,0,0,0,237,237,237,237,255,0,240,30,2,0,0,0,240,240,240,240,255,0,232,30,2,0,0,0,242,242,242,242,255,0,224,30,2,0,0,0,245,245,245,245,255,0,200,30,2,0,0,0,247,247,247,247,255,0,192,30,2,0,0,0,250,250,250,250,255,0,184,30,2,0,0,0,252,252,252,252,255,0,144,182,1,0,85,15,255,240,255,240,255,0,160,30,2,0,85,15,255,240,255,240,255,0,144,30,2,0,85,15,238,224,238,224,255,0,128,30,2,0,85,14,205,193,205,193,255,0,112,30,2,0,85,14,139,131,139,131,255,0,184,180,1,0,233,150,255,255,105,180,255,0,96,30,2,0,234,145,255,255,110,180,255,0,48,30,2,0,235,141,238,238,106,167,255,0,32,30,2,0,236,135,205,205,96,144,255,0,16,30,2,0,234,148,139,139,58,98,255,0,184,178,1,0,0,140,205,205,92,92,255,0,248,29,2,0,0,148,255,255,106,106,255,0,232,29,2,0,0,148,238,238,99,99,255,0,208,29,2,0,0,149,205,205,85,85,255,0,192,29,2,0,0,148,139,139,58,58,255,0,240,176,1,0,194,255,130,75,0,130,255,0,168,150,2,0,42,0,255,255,255,254,0,0,0,175,1,0,42,15,255,255,255,240,255,0,168,29,2,0,42,15,255,255,255,240,255,0,160,29,2,0,42,15,238,238,238,224,255,0,152,29,2,0,42,14,205,205,205,193,255,0,136,29,2,0,42,14,139,139,139,131,255,0,216,171,1,0,38,106,240,240,230,140,255,0,112,29,2,0,39,112,255,255,246,143,255,0,96,29,2,0])
.concat([39,112,238,238,230,133,255,0,88,29,2,0,39,111,205,205,198,115,255,0,80,29,2,0,39,111,139,139,134,78,255,0,56,169,1,0,170,20,250,230,230,250,255,0,32,167,1,0,240,15,255,255,240,245,255,0,48,29,2,0,240,15,255,255,240,245,255,0,32,29,2,0,239,15,238,238,224,229,255,0,8,29,2,0,240,14,205,205,193,197,255,0,184,28,2,0,239,14,139,139,131,134,255,0,160,165,1,0,64,255,252,124,252,0,255,0,136,163,1,0,38,49,255,255,250,205,255,0,168,28,2,0,38,49,255,255,250,205,255,0,152,28,2,0,37,50,238,238,233,191,255,0,136,28,2,0,38,49,205,205,201,165,255,0,112,28,2,0,39,49,139,139,137,112,255,0,176,161,1,0,137,63,230,173,216,230,255,0,80,28,2,0,138,64,255,191,239,255,255,0,200,27,2,0,138,64,238,178,223,238,255,0,184,27,2,0,138,63,205,154,192,205,255,0,160,27,2,0,137,64,139,104,131,139,255,0,232,159,1,0,0,119,240,240,128,128,255,0,56,158,1,0,127,31,255,224,255,255,255,0,128,27,2,0,127,31,255,224,255,255,255,0,48,27,2,0,127,31,238,209,238,238,255,0,0,27,2,0,127,31,205,180,205,205,255,0,232,26,2,0,127,31,139,122,139,139,255,0,216,26,2,0,35,115,238,238,221,130,255,0,120,26,2,0,35,116,255,255,236,139,255,0,104,26,2,0,35,115,238,238,220,130,255,0,56,26,2,0,35,115,205,205,190,112,255,0,240,25,2,0,35,115,139,139,129,76,255,0,112,156,1,0,42,40,250,250,250,210,255,0,80,154,1,0,0,0,211,211,211,211,255,0,56,149,1,0,0,0,211,211,211,211,255,0,104,147,1,0,248,73,255,255,182,193,255,0,88,25,2,0,249,81,255,255,174,185,255,0,56,25,2,0,248,81,238,238,162,173,255,0,32,25,2,0,249,80,205,205,140,149,255,0,8,25,2,0,249,80,139,139,95,101,255,0,136,145,1,0,12,132,255,255,160,122,255,0,240,24,2,0,12,132,255,255,160,122,255,0,192,24,2,0,11,132,238,238,149,114,255,0,176,24,2,0,12,133,205,205,129,98,255,0,144,24,2,0,12,133,139,139,87,66,255,0,40,144,1,0,125,209,178,32,178,170,255,0,192,142,1,0,143,117,250,135,206,250,255,0,128,24,2,0,143,79,255,176,226,255,255,0,96,24,2,0,143,79,238,164,211,238,255,0,80,24,2,0,142,79,205,141,182,205,255,0,32,24,2,0,143,78,139,96,123,139,255,0,16,24,2,0,175,143,255,132,112,255,255,0,104,141,1,0,148,56,153,119,136,153,255,0,240,139,1,0,148,56,153,119,136,153,255,0,88,138,1,0,151,52,222,176,196,222,255,0,240,23,2,0,151,53,255,202,225,255,255,0,224,23,2,0,151,53,238,188,210,238,255,0,208,23,2,0,151,53,205,162,181,205,255,0,184,23,2,0,150,53,139,110,123,139,255,0,8,136,1,0,42,31,255,255,255,224,255,0,144,23,2,0,42,31,255,255,255,224,255,0,128,23,2,0,42,31,238,238,238,209,255,0,112,23,2,0,42,31,205,205,205,180,255,0,96,23,2,0,42,31,139,139,139,122,255,0,32,132,1,0,85,192,205,50,205,50,255,0,136,130,1,0,21,20,250,250,240,230,255,0,16,129,1,0,212,255,255,255,0,255,255,0,48,23,2,0,212,255,255,255,0,255,255,0,24,23,2,0,212,255,238,238,0,238,255,0,8,23,2,0,212,255,205,205,0,205,255,0,232,22,2,0,212,255,139,139,0,139,255,0,208,65,2,0,239,185,176,176,48,96,255,0,224,22,2,0,228,203,255,255,52,179,255,0,216,22,2,0,228,203,238,238,48,167,255,0,192,22,2,0,228,204,205,205,41,144,255,0,184,22,2,0,228,203,139,139,28,98,255,0,32,126,1,0,113,128,205,102,205,170,255,0,160,124,1,0,170,255,205,0,0,205,255,0,200,122,1,0,204,152,211,186,85,211,255,0,160,22,2,0,203,153,255,224,102,255,255,0,88,22,2,0,203,153,238,209,95,238,255,0,72,22,2,0,203,153,205,180,82,205,255,0,56,22,2,0,203,154,139,122,55,139,255,0,104,121,1,0,183,124,219,147,112,219,255,0,8,22,2,0,183,125,255,171,130,255,255,0,248,21,2,0,183,125,238,159,121,238,255,0,200,21,2,0,183,125,205,137,104,205,255,0,184,21,2,0,183,124,139,93,71,139,255,0,160,119,1,0,103,169,179,60,179,113,255,0,168,117,1,0,176,143,238,123,104,238,255,0,240,115,1,0,111,255,250,0,250,154,255,0,56,114,1,0,125,167,209,72,209,204,255,0,0,113,1,0,228,228,199,199,21,133,255,0,248,111,1,0,170,198,112,25,25,112,255,0,240,110,1,0,106,9,255,245,255,250,255,0,208,172,2,0,4,30,255,255,228,225,255,0,88,21,2,0,4,30,255,255,228,225,255,0,72,21,2,0,4,30,238,238,213,210,255,0,48,21,2,0,3,29,205,205,183,181,255,0,8,21,2,0,5,29,139,139,125,123,255,0,168,171,2,0,26,73,255,255,228,181,255,0,96,170,2,0,25,81,255,255,222,173,255,0,224,20,2,0,25,81,255,255,222,173,255,0,208,20,2,0,25,82,238,238,207,161,255,0,176,20,2,0,25,82,205,205,179,139,255,0,160,20,2,0,25,82,139,139,121,94,255,0,136,54,2,0,170,255,128,0,0,128,255,0,192,127,1,0,170,255,128,0,0,128,255,0,136,218,1,0,42,0,255,255,255,254,0,0,88,166,2,0,27,23,253,253,245,230,255,0,16,163,2,0,56,192,142,107,142,35,255,0,128,20,2,0,56,193,255,192,255,62,255,0,112,20,2,0,56,192,238,179,238,58,255,0,72,20,2,0,56,192,205,154,205,50,255,0,8,20,2,0,56,192,139,105,139,34,255,0,168,161,2,0,27,255,255,255,165,0,255,0,240,19,2,0,27,255,255,255,165,0,255,0,232,19,2,0,27,255,238,238,154,0,255,0,88,19,2,0,27,255,205,205,133,0,255,0,80,19,2,0,27,255,139,139,90,0,255,0,112,160,2,0,11,255,255,255,69,0,255,0,64,19,2,0,11,255,255,255,69,0,255,0,48,19,2,0,11,255,238,238,64,0,255,0,32,19,2,0,11,255,205,205,55,0,255,0,208,18,2,0,11,255,139,139,37,0,255,0,32,159,2,0,214,123,218,218,112,214,255,0,144,18,2,0,214,124,255,255,131,250,255,0,136,18,2,0,214,124,238,238,122,233,255,0,80,18,2,0,214,124,205,205,105,201,255,0,72,18,2,0,213,124,139,139,71,137,255,0,176,157,2,0,38,72,238,238,232,170,255,0,128,156,2,0,85,100,251,152,251,152,255,0,240,17,2,0,85,101,255,154,255,154,255,0,224,17,2,0,85,100,238,144,238,144,255,0,168,17,2,0,85,100,205,124,205,124,255,0,152,17,2,0,85,100,139,84,139,84,255,0,48,155,2,0,127,67,238,175,238,238,255,0,128,17,2,0,127,68,255,187,255,255,255,0,104,17,2,0,127,68,238,174,238,238,255,0,72,17,2,0,127,68,205,150,205,205,255,0,8,17,2,0,127,67,139,102,139,139,255,0,176,153,2,0,241,124,219,219,112,147,255,0,240,16,2,0,241,125,255,255,130,171,255,0,224,16,2,0,241,125,238,238,121,159,255,0,184,16,2,0,241,125,205,205,104,137,255,0,168,16,2,0,241,124,139,139,71,93,255,0,240,151,2,0,26,41,255,255,239,213,255,0,8,150,2,0,20,70,255,255,218,185,255,0,144,16,2,0,20,70,255,255,218,185,255,0,128,16,2,0,19,69,238,238,203,173,255,0,64,16,2,0,19,69,205,205,175,149,255,0,48,16,2,0,20,69,139,139,119,101,255,0,120,148,2,0,20,176,205,205,133,63,255,0,64,147,2,0,247,63,255,255,192,203,255,0,240,15,2,0,245,73,255,255,181,197,255,0,216,15,2,0,245,73,238,238,169,184,255,0,208,15,2,0,245,74,205,205,145,158,255,0,200,15,2,0,245,73,139,139,99,108,255,0,8,146,2,0,212,70,221,221,160,221,255,0,176,15,2,0,212,68,255,255,187,255,255,0,128,15,2,0,212,68,238,238,174,238,255,0,120,15,2,0,212,68,205,205,150,205,255,0,88,15,2,0,212,67,139,139,102,139,255,0,176,144,2,0,132,59,230,176,224,230,255,0,216,38,2,0,196,221,240,160,32,240,255,0,208,14,2,0,191,207,255,155,48,255,255,0,200,14,2,0,192,207,238,145,44,238,255,0,136,14,2,0,192,207,205,125,38,205,255,0,120,14,2,0,192,207,139,85,26,139,255,0,128,32,2,0,0,255,255,255,0,0,255,0,48,14,2,0,0,255,255,255,0,0,255,0,40,14,2,0,0,255,238,238,0,0,255,0,16,14,2,0,0,255,205,205,0,0,255,0,8,14,2,0,0,255,139,139,0,0,255,0,184,140,2,0,0,61,188,188,143,143,255,0,176,13,2,0,0,62,255,255,193,193,255,0,160,13,2,0,0,62,238,238,180,180,255,0,144,13,2,0,0,62,205,205,155,155,255,0,112,13,2,0,0,62,139,139,105,105,255,0,136,139,2,0,159,181,225,65,105,225,255,0,40,13,2,0,159,183,255,72,118,255,255,0,24,13,2,0,159,183,238,67,110,238,255,0,8,13,2,0,159,182,205,58,95,205,255,0,248,12,2,0,159,183,139,39,64,139,255,0,176,137,2,0,17,220,139,139,69,19,255,0,8,136,2,0,4,138,250,250,128,114,255,0,144,12,2,0,9,150,255,255,140,105,255,0,128,12,2,0,9,150,238,238,130,98,255,0,112,12,2,0,9,150,205,205,112,84,255,0,104,12,2,0,9,150,139,139,76,57,255,0,208,134,2,0,19,154,244,244,164,96,255,0,112,133,2,0,103,170,139,46,139,87,255,0,64,12,2,0,103,171,255,84,255,159,255,0,48,12,2,0,103,171,238,78,238,148,255,0,16,12,2,0,103,171,205,67,205,128,255,0,0,12,2,0,103,170,139,46,139,87,255,0,184,131,2,0,17,16,255,255,245,238,255,0,240,11,2,0,17,16,255,255,245,238,255,0,216,11,2,0,18,17,238,238,229,222,255,0,176,11,2,0,18,17,205,205,197,191,255,0,112,11,2,0,18,16,139,139,134,130,255,0,160,130,2,0,13,183,160,160,82,45,255,0,104,11,2,0,13,184,255,255,130,71,255,0,96,11,2,0,13,184,238,238,121,66,255,0,72,11,2,0,13,184,205,205,104,57,255,0,64,11,2,0,13,185,139,139,71,38,255,0,248,127,2,0,139,108,235,135,206,235,255,0,48,11,2,0,144,120,255,135,206,255,255,0,24,11,2,0,144,120,238,126,192,238,255,0,8,11,2,0,144,120,205,108,166,205,255,0,224,10,2,0,145,119,139,74,112,139,255,0,104,126,2,0,175,143,205,106,90,205,255,0,200,10,2,0,175,144,255,131,111,255,255,0,184,10,2,0,175,144,238,122,103,238,255,0,152,10,2,0,175,144,205,105,89,205,255,0,120,10,2,0,175,144,139,71,60,139,255,0,160,124,2,0,148,56,144,112,128,144,255,0,24,10,2,0,149,56,255,198,226,255,255,0,248,9,2,0,149,56,238,185,211,238,255,0,232,9,2,0,148,57,205,159,182,205,255,0,168,9,2,0,149,56,139,108,123,139,255,0,208,122,2,0,148,56,144,112,128,144,255,0,80,121,2,0,0,5,255,255,250,250,255,0,120,9,2,0,0,5,255,255,250,250,255,0,24,9,2,0,0,5,238,238,233,233,255,0,16,9,2,0,0,4,205,205,201,201,255,0,8,9,2,0,0,3,139,139,137,137,255,0,8,120,2,0,106,255,255,0,255,127,255,0,176,8,2,0,106,255,255,0,255,127,255,0,160,8,2,0,106,255,238,0,238,118,255,0,72,8,2,0,106,255,205,0,205,102,255,0,32,8,2,0,106,255,139,0,139,69,255,0,176,118,2,0,146,155,180,70,130,180,255,0,200,7,2,0,146,156,255,99,184,255,255,0,112,7,2,0,146,156,238,92,172,238,255,0,96,7,2,0,146,156,205,79,148,205,255,0,80,7,2,0,147,155,139,54,100,139,255,0,160,116,2,0,24,84,210,210,180,140,255,0,48,7,2,0,20,176,255,255,165,79,255,0,0,7,2,0,20,176,238,238,154,73,255,0,224,6,2,0,20,176,205,205,133,63,255,0,216,6,2,0,20,176,139,139,90,43,255,0,240,113,2,0,212,29,216,216,191,216,255,0,176,6,2,0,212,30,255,255,225,255,255,0,136,6,2,0,212,30,238,238,210,238,255,0,120,6,2,0,212,29,205,205,181,205,255,0,104,6,2,0,212,29,139,139,123,139,255,0,200,112,2,0,6,184,255,255,99,71,255,0,80,6,2,0,6,184,255,255,99,71,255,0,72,6,2,0,6,184,238,238,92,66,255,0,48,6,2,0,6,184,205,205,79,57,255,0,40,6,2,0,6,185,139,139,54,38,255,0,112,123,1,0,42,0,255,255,255,254,0,0,104,111,2,0,123,182,224,64,224,208,255,0,176,5,2,0,129,255,255,0,245,255,255,0,152,5,2,0,129,255,238,0,229,238,255,0,136,5,2,0,129,255,205,0,197,205,255,0,120,5,2,0,129,255,139,0,134,139,255,0,144,109,2,0,212,115,238,238,130,238,255,0,0,142,2,0,227,215,208,208,32,144,255,0,72,5,2,0,235,193,255,255,62,150,255,0,56,5,2,0,235,192,238,238,58,140,255,0,40,5,2,0,235,192,205,205,50,120,255,0,24,5,2,0,235,192,139,139,34,82,255,0,96,108,2,0,27,68,245,245,222,179,255,0,0,5,2,0,27,69,255,255,231,186,255,0,248,4,2,0,27,68,238,238,216,174,255,0,240,4,2,0,27,68,205,205,186,150,255,0,224,4,2,0,27,67,139,139,126,102,255,0,192,6,2,0,0,0,255,255,255,255,255,0,96,105,2,0,0,0,245,245,245,245,255,0,104,2,2,0,42,255,255,255,255,0,255,0,216,4,2,0,42,255,255,255,255,0,255,0,208,4,2,0,42,255,238,238,238,0,255,0,184,4,2,0,42,255,205,205,205,0,255,0,176,4,2,0,42,255,139,139,139,0,255,0,240,102,2,0,56,192,205,154,205,50,255,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,81,160,79,228,73,210,14,64,180,200,118,190,159,58,53,192,58,34,223,165,212,37,213,191,243,130,62,71,154,46,138,63,159,229,121,112,119,214,249,191,126,253,16,27,44,156,230,63,150,236,216,8,196,235,204,63,205,206,162,119,42,224,208,63,176,227,191,64,16,32,237,191,173,161,212,94,68,219,216,63,59,161,124,230,81,150,118,63,211,110,112,249,122,132,123,63,129,204,206,162,119,42,228,191,209,173,215,244,160,160,200,63,106,223,55,25,176,63,132,63,190,202,144,25,94,255,132,63,28,150,6,126,84,195,196,191,165,73,41,232,246,226,35,64,169,217,3,173,192,144,193,63,8,196,144,65,147,105,137,63,250,68,158,36,93,51,208,191,1,240,153,54,45,194,94,63,13,156,125,47,207,148,151,63,137,181,248,20,0,227,137,63,229,169,88,70,52,203,177,191,143,0,201,207,161,103,166,191,92,181,198,251,204,180,136,63,77,164,143,84,58,179,144,63,230,199,4,161,97,214,160,191,199,105,103,28,19,247,130,191,42,127,107,229,45,112,92,191,228,87,98,84,8,154,117,63,209,241,135,85,114,4,183,63,149,212,9,104,34,60,51,192,100,35,16,175,235,119,16,192,167,33,170,240,103,120,199,63,218,255,0,107,213,174,193,63,78,40,68,192,33,84,247,191,170,72,133,177,133,32,245,63,157,104,87,33,229,39,246,63,77,46,198,192,58,142,205,63,89,107,40,181,23,209,220,191,3,63,170,97,191,39,204,63,166,71,83,61,153,127,218,63,182,129,59,80,167,60,174,63,81,76,222,0,51,223,185,191,245,118,149,255,218,11,166,63,212,165,53,188,15,246,148,63,31,173,32,188,44,220,144,63,40,44,241,128,178,201,35,64,35,90,225,76,2,138,183,63,72,163,101,81,150,41,127,63,187,180,134,247,193,158,147,63,23,168,123,83,71,125,160,191,33,43,174,224,109,148,139,63,51,115,220,132,214,30,181,191,160,120,132,137,245,252,143,63,105,53,36,238,177,244,145,191,184,205,51,122,94,191,106,63,146,62,173,162,63,52,205,191,126,176,231,198,79,62,152,191,7,35,155,80,45,199,164,63,62,24,194,123,88,185,145,191,45,124,125,173,75,141,198,63,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,76,1,0,0,162,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
.concat([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
.concat([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
.concat([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,244,108,86,125,174,182,214,63,181,21,251,203,238,201,225,63,181,21,251,203,238,201,225,63,196,66,173,105,222,113,236,63,16,122,54,171,62,87,229,63,245,219,215,129,115,70,204,63,88,168,53,205,59,78,213,63,88,168,53,205,59,78,213,63,136,133,90,211,188,227,216,63,1,77,132,13,79,175,226,63,211,188,227,20,29,201,209,63,88,168,53,205,59,78,213,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,181,21,251,203,238,201,225,63,181,21,251,203,238,201,225,63,181,21,251,203,238,201,225,63,181,21,251,203,238,201,225,63,181,21,251,203,238,201,225,63,181,21,251,203,238,201,225,63,181,21,251,203,238,201,225,63,181,21,251,203,238,201,225,63,181,21,251,203,238,201,225,63,181,21,251,203,238,201,225,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,1,77,132,13,79,175,226,63,1,77,132,13,79,175,226,63,1,77,132,13,79,175,226,63,181,21,251,203,238,201,225,63,204,93,75,200,7,61,240,63,16,122,54,171,62,87,229,63,16,122,54,171,62,87,229,63,210,111,95,7,206,25,231,63,210,111,95,7,206,25,231,63,16,122,54,171,62,87,229,63,120,11,36,40,126,140,227,63,106,222,113,138,142,228,232,63,210,111,95,7,206,25,231,63,211,188,227,20,29,201,209,63,0,0,0,0,0,0,224,63,16,122,54,171,62,87,229,63,181,21,251,203,238,201,225,63,44,212,154,230,29,167,234,63,210,111,95,7,206,25,231,63,106,222,113,138,142,228,232,63,16,122,54,171,62,87,229,63,106,222,113,138,142,228,232,63,210,111,95,7,206,25,231,63,16,122,54,171,62,87,229,63,120,11,36,40,126,140,227,63,210,111,95,7,206,25,231,63,16,122,54,171,62,87,229,63,134,56,214,197,109,52,238,63,16,122,54,171,62,87,229,63,16,122,54,171,62,87,229,63,120,11,36,40,126,140,227,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,166,10,70,37,117,2,222,63,181,21,251,203,238,201,225,63,72,191,125,29,56,103,204,63,181,21,251,203,238,201,225,63,181,21,251,203,238,201,225,63,0,0,0,0,0,0,224,63,181,21,251,203,238,201,225,63,181,21,251,203,238,201,225,63,211,188,227,20,29,201,209,63,181,21,251,203,238,201,225,63,181,21,251,203,238,201,225,63,72,191,125,29,56,103,204,63,72,191,125,29,56,103,204,63,0,0,0,0,0,0,224,63,72,191,125,29,56,103,204,63,44,212,154,230,29,167,234,63,181,21,251,203,238,201,225,63,181,21,251,203,238,201,225,63,181,21,251,203,238,201,225,63,181,21,251,203,238,201,225,63,88,168,53,205,59,78,213,63,0,0,0,0,0,0,224,63,211,188,227,20,29,201,209,63,181,21,251,203,238,201,225,63,0,0,0,0,0,0,224,63,210,111,95,7,206,25,231,63,0,0,0,0,0,0,224,63,0,0,0,0,0,0,224,63,0,0,0,0,0,0,224,63,2,154,8,27,158,94,213,63,224,190,14,156,51,162,208,63,2,154,8,27,158,94,213,63,1,77,132,13,79,175,226,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,88,168,53,205,59,78,213,63,181,21,251,203,238,201,225,63,181,21,251,203,238,201,225,63,62,232,217,172,250,92,197,63,181,21,251,203,238,201,225,63,181,21,251,203,238,201,225,63,181,21,251,203,238,201,225,63,181,21,251,203,238,201,225,63,29,56,103,68,105,111,200,63,88,168,53,205,59,78,213,63,181,21,251,203,238,201,225,63,88,168,53,205,59,78,213,63,88,168,53,205,59,78,213,63,0,0,0,0,0,0,224,63,0,0,0,0,0,0,224,63,211,188,227,20,29,201,209,63,181,21,251,203,238,201,225,63,181,21,251,203,238,201,225,63,181,21,251,203,238,201,225,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,231,29,167,232,72,46,225,63,162,180,55,248,194,100,214,63,72,191,125,29,56,103,204,63,88,168,53,205,59,78,213,63,88,168,53,205,59,78,213,63,181,21,251,203,238,201,225,63,0,0,0,0,0,0,240,63,0,0,0,0,0,0,240,63,211,188,227,20,29,201,209,63,120,11,36,40,126,140,227,63,211,188,227,20,29,201,209,63,88,168,53,205,59,78,213,63,88,168,53,205,59,78,213,63,88,168,53,205,59,78,213,63,88,168,53,205,59,78,213,63,88,168,53,205,59,78,213,63,88,168,53,205,59,78,213,63,88,168,53,205,59,78,213,63,88,168,53,205,59,78,213,63,211,188,227,20,29,201,209,63,88,168,53,205,59,78,213,63,88,168,53,205,59,78,213,63,211,188,227,20,29,201,209,63,88,168,53,205,59,78,213,63,88,168,53,205,59,78,213,63,88,168,53,205,59,78,213,63,0,0,0,0,0,0,240,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,0,0,0,0,0,0,240,63,211,188,227,20,29,201,209,63,234,149,178,12,113,172,215,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,181,21,251,203,238,201,225,63,106,222,113,138,142,228,232,63,0,0,0,0,0,0,240,63,152,221,147,135,133,90,215,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,196,66,173,105,222,113,236,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,72,191,125,29,56,103,204,63,120,11,36,40,126,140,227,63,134,56,214,197,109,52,238,63,120,11,36,40,126,140,227,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,211,188,227,20,29,201,209,63,3,0,0,0,200,60,0,0,3,0,0,0,160,60,0,0,3,0,0,0,16,60,0,0,3,0,0,0,112,59,0,0,3,0,0,0,72,59,0,0,3,0,0,0,32,59,0,0,3,0,0,0,248,58,0,0,3,0,0,0,232,59,0,0,3,0,0,0,192,59,0,0,0,0,0,0,0,53,0,0,0,0,0,0,216,52,0,0,0,0,0,0,176,52,0,0,0,0,0,0,56,52,0,0,0,0,0,0,16,52,0,0,0,0,0,0,232,51,0,0,0,0,0,0,192,51,0,0,0,0,0,0,136,52,0,0,0,0,0,0,96,52,0,0,4,0,0,0,168,53,0,0,0,0,0,0,0,0,0,0,1,0,0,0,128,57,0,0,0,0,0,0,0,0,0,0,1,0,0,0,96,58,0,0,0,0,0,0,0,0,0,0,240,205,1,0,32,204,1,0,216,196,1,0,208,200,1,0,128,199,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,96,65,2,0,0,0,0,0,136,218,1,0,1,0,0,0,16,27,2,0,7,0,0,0,152,18,2,0,3,0,0,0,160,235,1,0,5,0,0,0,104,3,2,0,15,0,0,0,216,0,2,0,8,0,0,0,216,0,2,0,16,0,0,0,168,254,1,0,4,0,0,0,168,254,1,0,17,0,0,0,200,252,1,0,5,0,0,0,200,252,1,0,2,0,0,0,224,250,1,0,6,0,0,0,24,248,1,0,4,0,0,0,240,245,1,0,7,0,0,0,192,243,1,0,7,0,0,0,168,240,1,0,5,0,0,0,72,237,1,0,8,0,0,0,88,234,1,0,8,0,0,0,72,237,1,0,9,0,0,0,40,232,1,0,7,0,0,0,240,229,1,0,10,0,0,0,144,227,1,0,7,0,0,0,0,225,1,0,11,0,0,0,176,222,1,0,6,0,0,0,200,220,1,0,12,0,0,0,32,219,1,0,9,0,0,0,200,220,1,0,13,0,0,0,16,217,1,0,8,0,0,0,240,212,1,0,14,0,0,0,136,210,1,0,8,0,0,0,48,208,1,0,18,0,0,0,72,206,1,0,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,108,110,114,0,0,0,0,0,0,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,23,17,2,2,2,2,2,2,2,2,2,2,2,2,2,18,16,2,19,2,2,22,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,20,2,21,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,14,2,15,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,3,4,5,6,7,8,9,10,11,12,13,0,0,0,0,0,0,0,0,0,0,0,34,12,13,14,35,15,9,16,17,10,16,17,202,16,17,45,68,69,252,1,6,246,15,7,246,36,2,16,17,47,48,54,76,77,40,38,59,60,42,54,49,57,61,63,47,58,64,217,65,48,62,37,55,67,53,74,43,56,75,72,0,0,0,0,0,2,2,1,0,3,3,1,0,1,0,1,1,1,0,2,1,1,0,2,2,3,1,1,0,4,0,1,3,1,3,5,3,1,1,1,1,2,0,1,0,4,2,0,2,1,1,3,2,1,0,3,2,1,0,1,1,0,1,1,1,3,0,0,0,24,25,25,25,26,27,28,28,29,29,30,30,31,31,32,32,33,33,34,34,35,36,36,38,37,37,39,39,40,40,40,41,41,42,42,42,43,43,44,44,45,46,46,47,48,48,49,50,51,53,52,54,54,54,55,55,55,56,56,57,57,0,0,241,241,255,241,241,241,241,241,241,31,32,241,0,242,241,241,12,241,241,241,8,13,241,241,241,249,241,241,241,241,241,241,245,241,0,0,0,0,0,0,18,241,241,20,9,3,241,254,241,241,241,1,241,241,241,1,241,241,10,254,241,19,25,21,241,19,1,241,241,241,241,11,17,241,241,241,241,241,241,241,241,241,1,241,241,22,9,1,1,29,15,23,241,241,26,23,27,241,241,28,241,241,25,241,1,241,251,241,241,1,241,16,241,241,30,241,241,241,241,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,255,3,8,4,33,5,11,18,19,39,20,21,22,41,50,23,24,25,26,44,51,52,66,70,71,27,73,28,29,46,30,78,31,32,0,0,0,0,0,0,0,3,9,0,0,0,1,14,2,11,12,8,34,35,36,53,58,60,0,13,16,18,26,22,27,18,38,49,33,23,50,29,59,6,7,52,5,15,17,20,24,40,0,19,40,0,0,0,0,0,54,21,39,28,29,0,32,37,51,30,47,61,26,43,0,25,0,31,41,0,42,57,45,46,0,48,55,56,44,0,11,3,4,5,15,7,3,12,13,6,12,13,14,12,13,26,21,22,0,1,0,3,7,14,6,15,8,12,13,18,19,42,16,17,9,16,47,48,17,50,23,19,13,20,18,46,18,20,62,19,50,19,44,64,42,66,25,44,69,66,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10,0,11,0,12,0,13,0,14,0,10,0,15,0,16,0,17,0,18,0,19,0,10,0,20,0,21,0,21,0,21,0,22,0,23,0,21,0,24,0,21,0,21,0,25,0,21,0,21,0,21,0,26,0,21,0,21,0,10,0,21,0,21,0,21,0,22,0,23,0,24,0,21,0,21,0,25,0,21,0,21,0,21,0,26,0,21,0,21,0,12,0,12,0,35,0,29,0,29,0,31,0,32,0,31,0,32,0,35,0,36,0,37,0,44,0,49,0,46,0,45,0,41,0,36,0,37,0,39,0,40,0,50,0,41,0,51,0,42,0,52,0,53,0,54,0,58,0,49,0,48,0,59,0,33,0,67,0,33,0,61,0,62,0,68,0,50,0,51,0,69,0,52,0,53,0,54,0,46,0,89,0,41,0,43,0,89,0,67,0,66,0,70,0,72,0,68,0,71,0,89,0,58,0,69,0,89,0,59,0,73,0,89,0,63,0,66,0,74,0,43,0,75,0,76,0,70,0,72,0,71,0,77,0,78,0,79,0,43,0,80,0,73,0,81,0,89,0,82,0,83,0,74,0,75,0,84,0,76,0,85,0,27,0,77,0,78,0,86,0,79,0,80,0,87,0,88,0,81,0,82,0,83,0,89,0,89,0,84,0,89,0,89,0,85,0,89,0,89,0,86,0,89,0,89,0,87,0,88,0,28,0,28,0,28,0,28,0,28,0,28,0,28,0,30,0,30,0,30,0,30,0,30,0,30,0,30,0,34,0,34,0,34,0,34,0,34,0,34,0,34,0,38,0,89,0,38,0,38,0,38,0,38,0,38,0,47,0,47,0,55,0,89,0,55,0,55,0,55,0,55,0,55,0,56,0,89,0,56,0,89,0,56,0,56,0,56,0,57,0,89,0,57,0,57,0,57,0,57,0,57,0,60,0,60,0,89,0,60,0,60,0,60,0,60,0,64,0,89,0,64,0,64,0,64,0,64,0,65,0,89,0,65,0,65,0,65,0,65,0,65,0,9,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,2,0,0,0,3,0,0,0,1,0,0,0,4,0,0,0,1,0,0,0,5,0,0,0,1,0,0,0,6,0,0,0,7,0,0,0,7,0,0,0,1,0,0,0,6,0,0,0,6,0,0,0,6,0,0,0,6,0,0,0,6,0,0,0,6,0,0,0,6,0,0,0,6,0,0,0,6,0,0,0,6,0,0,0,6,0,0,0,6,0,0,0,6,0,0,0,6,0,0,0,6,0,0,0,6,0,0,0,3,0,0,0,6,0,0,0,6,0,0,0,6,0,0,0,6,0,0,0,6,0,0,0,6,0,0,0,6,0,0,0,6,0,0,0,6,0,0,0,6,0,0,0,6,0,0,0,6,0,0,0,6,0,0,0,6,0,0,0,6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,2,0,0,0,3,0,0,0,1,0,0,0,1,0,0,0,2,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,2,0,0,0,1,0,0,0,4,0,0,0,5,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,6,0,0,0,1,0,0,0,1,0,0,0,7,0,0,0,8,0,0,0,9,0,0,0,10,0,0,0,10,0,0,0,10,0,0,0,10,0,0,0,10,0,0,0,10,0,0,0,10,0,0,0,10,0,0,0,10,0,0,0,10,0,0,0,1,0,0,0,1,0,0,0,11,0,0,0,1,0,0,0,12,0,0,0,1,0,0,0,13,0,0,0,14,0,0,0,15,0,0,0,16,0,0,0,17,0,0,0,18,0,0,0,19,0,0,0,20,0,0,0,21,0,0,0,22,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,23,0,0,0,24,0,0,0,25,0,0,0,19,0,0,0,26,0,0,0,27,0,0,0,28,0,0,0,29,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,1,0,0,0,30,0,0,0,1,0,0,0,1,0,0,0,19,0,0,0,1,0,0,0,31,0,0,0,32,0,0,0,33,0,0,0,34,0,0,0,35,0,0,0,19,0,0,0,36,0,0,0,37,0,0,0,38,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,39,0,0,0,40,0,0,0,41,0,0,0,19,0,0,0,42,0,0,0,43,0,0,0,44,0,0,0,45,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0])
.concat([19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,19,0,0,0,0,0,89,0,1,0,90,0,90,0,91,0,91,0,92,0,92,0,89,0,89,0,89,0,89,0,89,0,93,0,89,0,89,0,89,0,94,0,89,0,89,0,95,0,95,0,95,0,95,0,95,0,95,0,96,0,97,0,98,0,99,0,99,0,89,0,89,0,100,0,89,0,89,0,89,0,93,0,89,0,89,0,94,0,89,0,94,0,89,0,101,0,94,0,89,0,95,0,95,0,95,0,95,0,95,0,95,0,95,0,96,0,97,0,98,0,98,0,89,0,99,0,89,0,89,0,89,0,100,0,101,0,94,0,95,0,95,0,95,0,95,0,95,0,95,0,95,0,95,0,95,0,95,0,95,0,95,0,95,0,95,0,95,0,95,0,95,0,95,0,95,0,95,0,95,0,95,0,0,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,0,0,0,0,0,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,3,0,4,0,7,0,3,0,4,0,5,0,5,0,6,0,6,0,8,0,7,0,7,0,17,0,22,0,18,0,17,0,18,0,8,0,8,0,15,0,15,0,23,0,15,0,24,0,15,0,25,0,26,0,26,0,29,0,22,0,95,0,29,0,5,0,49,0,6,0,33,0,33,0,50,0,23,0,24,0,51,0,25,0,26,0,26,0,41,0,43,0,41,0,43,0,46,0,49,0,46,0,52,0,54,0,50,0,53,0,57,0,58,0,51,0,57,0,58,0,67,0,66,0,33,0,66,0,68,0,40,0,69,0,70,0,52,0,54,0,53,0,71,0,72,0,73,0,16,0,75,0,67,0,77,0,9,0,78,0,79,0,68,0,69,0,81,0,70,0,82,0,2,0,71,0,72,0,83,0,73,0,75,0,85,0,87,0,77,0,78,0,79,0,0,0,0,0,81,0,0,0,0,0,82,0,0,0,0,0,83,0,0,0,0,0,85,0,87,0,90,0,90,0,90,0,90,0,90,0,90,0,90,0,91,0,91,0,91,0,91,0,91,0,91,0,91,0,92,0,92,0,92,0,92,0,92,0,92,0,92,0,93,0,0,0,93,0,93,0,93,0,93,0,93,0,94,0,94,0,96,0,0,0,96,0,96,0,96,0,96,0,96,0,97,0,0,0,97,0,0,0,97,0,97,0,97,0,98,0,0,0,98,0,98,0,98,0,98,0,98,0,99,0,99,0,0,0,99,0,99,0,99,0,99,0,100,0,0,0,100,0,100,0,100,0,100,0,101,0,0,0,101,0,101,0,101,0,101,0,101,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,89,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,127,0,43,0,44,0,48,0,50,0,45,0,52,0,124,0,227,0,227,0,227,0,227,0,0,0,58,0,110,0,52,0,52,0,227,0,227,0,0,0,37,0,50,0,43,0,47,0,44,0,0,0,0,0,68,0,0,0,0,0,227,0,78,0,0,0,227,0,227,0,227,0,0,0,227,0,101,0,82,0,227,0,83,0,227,0,0,0,86,0,227,0,0,0,59,0,63,0,72,0,80,0,74,0,83,0,0,0,0,0,95,0,96,0,227,0,0,0,227,0,227,0,227,0,0,0,0,0,99,0,80,0,92,0,87,0,95,0,95,0,98,0,105,0,0,0,100,0,0,0,107,0,99,0,101,0,0,0,101,0,117,0,114,0,0,0,113,0,0,0,118,0,0,0,227,0,155,0,162,0,169,0,176,0,179,0,70,0,185,0,192,0,199,0,206,0,213,0,219,0,0,0,0,0,0,0,0,0,0,0,4,0,4,0,26,0,26,0,31,0,31,0,34,0,32,0,10,0,2,0,21,0,9,0,32,0,32,0,32,0,20,0,27,0,1,0,19,0,19,0,19,0,19,0,19,0,19,0,8,0,4,0,5,0,26,0,2,0,22,0,26,0,31,0,30,0,29,0,28,0,9,0,18,0,0,0,20,0,17,0,20,0,3,0,7,0,20,0,20,0,19,0,19,0,19,0,19,0,19,0,19,0,19,0,8,0,4,0,5,0,5,0,6,0,26,0,25,0,23,0,24,0,31,0,7,0,20,0,19,0,19,0,19,0,19,0,19,0,19,0,19,0,12,0,19,0,11,0,19,0,19,0,19,0,13,0,19,0,19,0,19,0,15,0,19,0,14,0,19,0,16,0,0,0,0,0,0,0,47,114,100,98,117,52,47,50,0,0,0,0,0,0,0,0,47,114,100,98,117,52,47,49,0,0,0,0,0,0,0,0,47,114,100,98,117,51,47,51,0,0,0,0,0,0,0,0,74,117,108,0,0,0,0,0,47,114,100,98,117,51,47,50,0,0,0,0,0,0,0,0,32,32,32,60,83,67,82,73,80,84,32,76,65,78,71,85,65,71,69,61,39,74,97,118,97,115,99,114,105,112,116,39,62,10,0,0,0,0,0,0,47,114,100,98,117,51,47,49,0,0,0,0,0,0,0,0,109,105,110,116,99,114,101,97,109,0,0,0,0,0,0,0,47,115,101,116,95,115,99,97,108,101,32,123,0,0,0,0,47,114,100,98,117,49,49,47,57,0,0,0,0,0,0,0,90,101,116,97,0,0,0,0,47,114,100,98,117,49,49,47,56,0,0,0,0,0,0,0,100,101,99,111,114,97,116,101,0,0,0,0,0,0,0,0,101,100,103,101,116,97,114,103,101,116,0,0,0,0,0,0,47,114,100,98,117,49,49,47,55,0,0,0,0,0,0,0,47,114,100,98,117,49,49,47,54,0,0,0,0,0,0,0,47,114,100,98,117,49,49,47,53,0,0,0,0,0,0,0,47,98,114,98,103,51,47,50,0,0,0,0,0,0,0,0,47,114,100,98,117,49,49,47,52,0,0,0,0,0,0,0,47,114,100,98,117,49,49,47,51,0,0,0,0,0,0,0,47,114,100,98,117,49,49,47,50,0,0,0,0,0,0,0,74,117,110,0,0,0,0,0,47,114,100,98,117,49,49,47,49,49,0,0,0,0,0,0,60,33,45,45,32,80,97,103,101,115,58,32,37,100,32,45,45,62,10,0,0,0,0,0,109,105,100,110,105,103,104,116,98,108,117,101,0,0,0,0,47,73,110,118,83,99,97,108,101,70,97,99,116,111,114,32,49,46,48,32,100,101,102,0,47,114,100,98,117,49,49,47,49,48,0,0,0,0,0,0,47,114,100,98,117,49,49,47,49,0,0,0,0,0,0,0,89,117,109,108,0,0,0,0,47,114,100,98,117,49,48,47,57,0,0,0,0,0,0,0,109,105,110,108,101,110,0,0,116,97,114,103,101,116,0,0,47,114,100,98,117,49,48,47,56,0,0,0,0,0,0,0,47,114,100,98,117,49,48,47,55,0,0,0,0,0,0,0,47,114,100,98,117,49,48,47,54,0,0,0,0,0,0,0,47,98,114,98,103,51,47,49,0,0,0,0,0,0,0,0,47,114,100,98,117,49,48,47,53,0,0,0,0,0,0,0,47,114,100,98,117,49,48,47,52,0,0,0,0,0,0,0,47,114,100,98,117,49,48,47,51,0,0,0,0,0,0,0,65,112,114,0,0,0,0,0,47,114,100,98,117,49,48,47,50,0,0,0,0,0,0,0,60,47,84,73,84,76,69,62,0,0,0,0,0,0,0,0,109,101,100,105,117,109,118,105,111,108,101,116,114,101,100,0,47,99,111,111,114,100,102,111,110,116,32,99,111,111,114,100,45,102,111,110,116,45,102,97,109,105,108,121,32,102,105,110,100,102,111,110,116,32,56,32,115,99,97,108,101,102,111,110,116,32,100,101,102,0,0,0,47,114,100,98,117,49,48,47,49,48,0,0,0,0,0,0,47,114,100,98,117,49,48,47,49,0,0,0,0,0,0,0,114,0,0,0,0,0,0,0,89,97,99,117,116,101,0,0,109,97,114,103,105,110,0,0,47,112,117,114,112,108,101,115,57,47,57,0,0,0,0,0,108,97,98,101,108,97,110,103,108,101,0,0,0,0,0,0,104,101,97,100,85,82,76,0,47,112,117,114,112,108,101,115,57,47,56,0,0,0,0,0,47,112,117,114,112,108,101,115,57,47,55,0,0,0,0,0,47,112,117,114,112,108,101,115,57,47,54,0,0,0,0,0,47,98,114,98,103,49,49,47,57,0,0,0,0,0,0,0,47,112,117,114,112,108,101,115,57,47,53,0,0,0,0,0,47,112,117,114,112,108,101,115,57,47,52,0,0,0,0,0,47,112,117,114,112,108,101,115,57,47,51,0,0,0,0,0,77,97,114,0,0,0,0,0,47,112,117,114,112,108,101,115,57,47,50,0,0,0,0,0,60,84,73,84,76,69,62,0,109,101,100,105,117,109,116,117,114,113,117,111,105,115,101,0,47,100,101,102,97,117,108,116,45,102,111,110,116,45,102,97,109,105,108,121,32,47,84,105,109,101,115,45,82,111,109,97,110,32,100,101,102,0,0,0,47,112,117,114,112,108,101,115,57,47,49,0,0,0,0,0,47,97,99,99,101,110,116,53,47,49,0,0,0,0,0,0,47,112,117,114,112,108,101,115,56,47,56,0,0,0,0,0,88,105,0,0,0,0,0,0,92,78,0,0,0,0,0,0,47,112,117,114,112,108,101,115,56,47,55,0,0,0,0,0,108,97,98,101,108,100,105,115,116,97,110,99,101,0,0,0,104,101,97,100,104,114,101,102,0,0,0,0,0,0,0,0,112,101,110,0,0,0,0,0,47,112,117,114,112,108,101,115,56,47,54,0,0,0,0,0,47,112,117,114,112,108,101,115,56,47,53,0,0,0,0,0,47,112,117,114,112,108,101,115,56,47,52,0,0,0,0,0,47,98,114,98,103,49,49,47,56,0,0,0,0,0,0,0,117,110,115,117,112,112,111,114,116,101,100,32,108,111,99,97,108,101,32,102,111,114,32,115,116,97,110,100,97,114,100,32,105,110,112,117,116,0,0,0,108,111,110,103,0,0,0,0,47,112,117,114,112,108,101,115,56,47,51,0,0,0,0,0,47,112,117,114,112,108,101,115,56,47,50,0,0,0,0,0,118,109,108,58,118,109,108,0,47,112,117,114,112,108,101,115,56,47,49,0,0,0,0,0,70,101,98,0,0,0,0,0,47,112,117,114,112,108,101,115,55,47,55,0,0,0,0,0,60,77,69,84,65,32,104,116,116,112,45,101,113,117,105,118,61,34,67,111,110,116,101,110,116,45,84,121,112,101,34,32,99,111,110,116,101,110,116,61,34,116,101,120,116,47,104,116,109,108,59,32,99,104,97,114,115,101,116,61,85,84,70,45,56,34,62,10,0,0,0,0,109,101,100,105,117,109,115,112,114,105,110,103,103,114,101,101,110,0,0,0,0,0,0,0,47,99,111,111,114,100,45,102,111,110,116,45,102,97,109,105,108,121,32,47,84,105,109,101,115,45,82,111,109,97,110,32,100,101,102,0,0,0,0,0,47,112,117,114,112,108,101,115,55,47,54,0,0,0,0,0,47,112,117,114,112,108,101,115,55,47,53,0,0,0,0,0,85,117,109,108,0,0,0,0,98,97,100,32,108,97,98,101,108,32,102,111,114,109,97,116,32,37,115,10,0,0,0,0,47,112,117,114,112,108,101,115,55,47,52,0,0,0,0,0,108,97,98,101,108,102,111,110,116,99,111,108,111,114,0,0,116,97,105,108,85,82,76,0,115,101,116,108,105,110,101,119,105,100,116,104,0,0,0,0,47,112,117,114,112,108,101,115,55,47,51,0,0,0,0,0,117,110,102,105,108,108,101,100,0,0,0,0,0,0,0,0,47,112,117,114,112,108,101,115,55,47,50,0,0,0,0,0,98,114,97,115,115,0,0,0,47,112,117,114,112,108,101,115,55,47,49,0,0,0,0,0,47,98,114,98,103,49,49,47,55,0,0,0,0,0,0,0,102,105,108,108,32,0,0,0,112,111,108,121,32,37,115,0,47,112,117,114,112,108,101,115,54,47,54,0,0,0,0,0,112,108,97,105,110,45,101,120,116,58,100,111,116,0,0,0,47,112,117,114,112,108,101,115,54,47,53,0,0,0,0,0,106,112,103,58,102,105,103,0,110,111,112,50,0,0,0,0,65,103,101,100,103,101,105,110,102,111,95,116,0,0,0,0,47,112,117,114,112,108,101,115,54,47,52,0,0,0,0,0,108,97,98,101,108,95,115,99,104,101,109,101,0,0,0,0,74,97,110,0,0,0,0,0,33,102,108,97,103,0,0,0,47,112,117,114,112,108,101,115,54,47,51,0,0,0,0,0,60,72,69,65,68,62,0,0,109,101,100,105,117,109,115,108,97,116,101,98,108,117,101,0,37,37,66,101,103,105,110,82,101,115,111,117,114,99,101,58,32,112,114,111,99,115,101,116,32,103,114,97,112,104,118,105,122,32,48,32,48,0,0,0,85,115,105,110,103,32,37,115,58,32,37,115,58,37,115,10,0,0,0,0,0,0,0,0,47,112,117,114,112,108,101,115,54,47,50,0,0,0,0,0,47,112,117,114,112,108,101,115,54,47,49,0,0,0,0,0,113,0,0,0,0,0,0,0,85,112,115,105,108,111,110,0,110,111,100,101,32,37,115,32,105,110,32,103,114,97,112,104,32,37,115,32,104,97,115,32,110,111,32,112,111,115,105,116,105,111,110,10,0,0,0,0,47,112,117,114,112,108,101,115,53,47,53,0,0,0,0,0,108,97,98,101,108,102,111,110,116,110,97,109,101,0,0,0,116,97,105,108,104,114,101,102,0,0,0,0,0,0,0,0,47,112,117,114,112,108,101,115,53,47,52,0,0,0,0,0,47,112,117,114,112,108,101,115,53,47,51,0,0,0,0,0,115,101,112,0,0,0,0,0,112,108,117,115,0,0,0,0,47,112,117,114,112,108,101,115,53,47,50,0,0,0,0,0,47,98,114,98,103,49,49,47,54,0,0,0,0,0,0,0,127,98,111,116,0,0,0,0,104,112,0,0,0,0,0,0,111,117,116,0,0,0,0,0,47,112,117,114,112,108,101,115,53,47,49,0,0,0,0,0,123,37,115,125,0,0,0,0,37,115,32,45,62,32,37,115,58,32,104,101,97,100,32,105,115,32,105,110,115,105,100,101,32,116,97,105,108,32,99,108,117,115,116,101,114,32,37,115,10,0,0,0,0,0,0,0,47,112,117,114,112,108,101,115,52,47,52,0,0,0,0,0,47,112,117,114,112,108,101,115,52,47,51,0,0,0,0,0,99,97,110,110,111,116,32,109,97,108,108,111,99,32,112,110,108,115,0,0,0,0,0,0,68,101,99,101,109,98,101,114,0,0,0,0,0,0,0,0,47,112,117,114,112,108,101,115,52,47,50,0,0,0,0,0,60,47,66,79,68,89,62,10,60,47,72,84,77,76,62,10,0,0,0,0,0,0,0,0,109,101,100,105,117,109,115,101,97,103,114,101,101,110,0,0,125,32,98,105,110,100,32,100,101,102,0,0,0,0,0,0,111,108,100,103,111,108,100,0,47,112,117,114,112,108,101,115,52,47,49,0,0,0,0,0,47,112,117,114,112,108,101,115,51,47,51,0,0,0,0,0,60,110,117,108,108,62,0,0,85,103,114,97,118,101,0,0,105,110,32,110,111,100,101,32,37,115,10,0,0,0,0,0,47,112,117,114,112,108,101,115,51,47,50,0,0,0,0,0,108,97,98,101,108,102,111,110,116,115,105,122,101,0,0,0,108,97,98,101,108,85,82,76,0,0,0,0,0,0,0,0,47,112,117,114,112,108,101,115,51,47,49,0,0,0,0,0,105,110,32,99,104,101,99,107,112,97,116,104,44,32,98,111,120,32,48,32,104,97,115,32,76,76,32,99,111,111,114,100,32,62,32,85,82,32,99,111,111,114,100,10,0,0,0,0,84,82,65,73,76,69,82,0,47,112,117,114,100,57,47,57,0,0,0,0,0,0,0,0,116,114,111,117,98,108,101,32,105,110,32,105,110,105,116,95,114,97,110,107,10,0,0,0,47,112,117,114,100,57,47,56,0,0,0,0,0,0,0,0,47,98,114,98,103,49,49,47,53,0,0,0,0,0,0,0,47,112,117,114,100,57,47,55,0,0,0,0,0,0,0,0,61,0,0,0,0,0,0,0,47,112,117,114,100,57,47,54,0,0,0,0,0,0,0,0,47,112,117,114,100,57,47,53,0,0,0,0,0,0,0,0,85,82,87,32,71,111,116,104,105,99,32,76,0,0,0,0,78,111,118,101,109,98,101,114,0,0,0,0,0,0,0,0,47,112,117,114,100,57,47,52,0,0,0,0,0,0,0,0,60,33,45,45,32,105,110,115,101,114,116,32,97,110,121,32,111,116,104,101,114,32,78,79,78,45,73,69,32,104,116,109,108,32,99,111,110,116,101,110,116,32,104,101,114,101,32,45,45,62,10,0,0,0,0,0,109,101,100,105,117,109,112,117,114,112,108,101,0,0,0,0,99,108,101,97,114,116,111,109,97,114,107,0,0,0,0,0,110,101,119,116,97,110,0,0,47,112,117,114,100,57,47,51,0,0,0,0,0,0,0,0,47,112,117,114,100,57,47,50,0,0,0,0,0,0,0,0,85,99,105,114,99,0,0,0,102,108,101,120,32,115,99,97,110,110,101,114,32,112,117,115,104,45,98,97,99,107,32,111,118,101,114,102,108,111,119,0,47,112,117,114,100,57,47,49,0,0,0,0,0,0,0,0,116,97,105,108,108,97,98,101,108,0,0,0,0,0,0,0,108,97,98,101,108,104,114,101,102,0,0,0,0,0,0,0,47,112,117,114,100,56,47,56,0,0,0,0,0,0,0,0,47,112,117,114,100,56,47,55,0,0,0,0,0,0,0,0,47,112,117,114,100,56,47,54,0,0,0,0,0,0,0,0,47,98,114,98,103,49,49,47,52,0,0,0,0,0,0,0,47,112,117,114,100,56,47,53,0,0,0,0,0,0,0,0,47,112,117,114,100,56,47,52,0,0,0,0,0,0,0,0,47,112,117,114,100,56,47,51,0,0,0,0,0,0,0,0,79,99,116,111,98,101,114,0,47,112,117,114,100,56,47,50,0,0,0,0,0,0,0,0,60,68,73,86,32,105,100,61,39,95,110,111,116,86,77,76,50,95,39,32,115,116,121,108,101,61,34,112,111,115,105,116,105,111,110,58,114,101,108,97,116,105,118,101,59,34,62,10,0,0,0,0,0,0,0,0,109,101,100,105,117,109,111,114,99,104,105,100,0,0,0,0,47,67,111,117,114,105,101,114,45,66,111,108,100,79,98,108,105,113,117,101,32,115,116,97,114,110,101,116,73,83,79,32,100,101,102,0,0,0,0,0,110,101,119,109,105,100,110,105,103,104,116,98,108,117,101,0,37,115,32,115,97,118,101,32,112,111,105,110,116,32,115,105,122,101,32,97,110,100,32,102,111,110,116,10,46,110,114,32,46,83,32,92,110,40,46,115,10,46,110,114,32,68,70,32,92,110,40,46,102,10,0,0,47,112,117,114,100,56,47,49,0,0,0,0,0,0,0,0,47,112,117,114,100,55,47,55,0,0,0,0,0,0,0,0,85,97,99,117,116,101,0,0,116,114,97,110,115,112,97,114,101,110,116,0,0,0,0,0,47,112,117,114,100,55,47,54,0,0,0,0,0,0,0,0,104,101,97,100,108,97,98,101,108,0,0,0,0,0,0,0,101,100,103,101,85,82,76,0,47,112,117,114,100,55,47,53,0,0,0,0,0,0,0,0,47,112,117,114,100,55,47,52,0,0,0,0,0,0,0,0,47,112,117,114,100,55,47,51,0,0,0,0,0,0,0,0,47,98,114,98,103,49,49,47,51,0,0,0,0,0,0,0,47,112,117,114,100,55,47,50,0,0,0,0,0,0,0,0,47,112,117,114,100,55,47,49,0,0,0,0,0,0,0,0,47,112,117,114,100,54,47,54,0,0,0,0,0,0,0,0,83,101,112,116,101,109,98,101,114,0,0,0,0,0,0,0,47,112,117,114,100,54,47,53,0,0,0,0,0,0,0,0,60,72,50,62,83,111,114,114,121,44,32,116,104,105,115,32,100,105,97,103,114,97,109,32,119,105,108,108,32,111,110,108,121,32,100,105,115,112,108,97,121,32,99,111,114,114,101,99,116,108,121,32,111,110,32,73,110,116,101,114,110,101,116,32,69,120,112,108,111,114,101,114,32,53,32,40,97,110,100,32,117,112,41,32,98,114,111,119,115,101,114,115,46,60,47,72,50,62,10,0,0,0,0,0,109,101,100,105,117,109,98,108,117,101,0,0,0,0,0,0,47,67,111,117,114,105,101,114,45,66,111,108,100,32,115,116,97,114,110,101,116,73,83,79,32,100,101,102,0,0,0,0,110,101,111,110,112,105,110,107,0,0,0,0,0,0,0,0,37,115,32,84,105,116,108,101,58,32,37,115,10,0,0,0,47,112,117,114,100,54,47,52,0,0,0,0,0,0,0,0,47,112,117,114,100,54,47,51,0,0,0,0,0,0,0,0,84,104,101,116,97,0,0,0,102,105,108,108,101,100,0,0,47,112,117,114,100,54,47,50,0,0,0,0,0,0,0,0,97,114,114,111,119,116,97,105,108,0,0,0,0,0,0,0,101,100,103,101,104,114,101,102,0,0,0,0,0,0,0,0,47,112,117,114,100,54,47,49,0,0,0,0,0,0,0,0,47,112,117,114,100,53,47,53,0,0,0,0,0,0,0,0,47,112,117,114,100,53,47,52,0,0,0,0,0,0,0,0,47,98,114,98,103,49,49,47,50,0,0,0,0,0,0,0,103,118,99,66,117,105,108,100,68,97,116,101,0,0,0,0,47,112,117,114,100,53,47,51,0,0,0,0,0,0,0,0,47,112,117,114,100,53,47,50,0,0,0,0,0,0,0,0,47,112,117,114,100,53,47,49,0,0,0,0,0,0,0,0,65,117,103,117,115,116,0,0,47,112,117,114,100,52,47,52,0,0,0,0,0,0,0,0,60,33,45,45,32,116,104,105,115,32,115,104,111,117,108,100,32,111,110,108,121,32,100,105,115,112,108,97,121,32,111,110,32,78,79,78,45,73,69,32,98,114,111,119,115,101,114,115,32,45,45,62,10,0,0,0,109,101,100,105,117,109,97,113,117,97,109,97,114,105,110,101,0,0,0,0,0,0,0,0,47,67,111,117,114,105,101,114,45,79,98,108,105,113,117,101,32,115,116,97,114,110,101,116,73,83,79,32,100,101,102,0,110,101,111,110,98,108,117,101,0,0,0,0,0,0,0,0,37,115,32,67,114,101,97,116,111,114,58,32,37,115,32,118,101,114,115,105,111,110,32,37,115,32,40,37,115,41,10,0,47,112,117,114,100,52,47,51,0,0,0,0,0,0,0,0,47,112,117,114,100,52,47,50,0,0,0,0,0,0,0,0,111,110,101,0,0,0,0,0,84,97,117,0,0,0,0,0,105,110,118,105,115,0,0,0,47,112,117,114,100,52,47,49,0,0,0,0,0,0,0,0,97,114,114,111,119,104,101,97,100,0,0,0,0,0,0,0,85,82,76,0,0,0,0,0,47,112,117,114,100,51,47,51,0,0,0,0,0,0,0,0,47,112,117,114,100,51,47,50,0,0,0,0,0,0,0,0,47,112,117,114,100,51,47,49,0,0,0,0,0,0,0,0,47,98,114,98,103,49,49,47,49,49,0,0,0,0,0,0,47,112,117,111,114,57,47,57,0,0,0,0,0,0,0,0,47,112,117,111,114,57,47,56,0,0,0,0,0,0,0,0,47,112,117,111,114,57,47,55,0,0,0,0,0,0,0,0,74,117,108,121,0,0,0,0,47,112,117,111,114,57,47,54,0,0,0,0,0,0,0,0,60,68,73,86,32,105,100,61,39,95,110,111,116,86,77,76,49,95,39,32,115,116,121,108,101,61,34,112,111,115,105,116,105,111,110,58,114,101,108,97,116,105,118,101,59,34,62,10,0,0,0,0,0,0,0,0,47,67,111,117,114,105,101,114,32,115,116,97,114,110,101,116,73,83,79,32,100,101,102,0,110,97,118,121,98,108,117,101,0,0,0,0,0,0,0,0,37,115,32,114,101,115,116,111,114,101,32,112,111,105,110,116,32,115,105,122,101,32,97,110,100,32,102,111,110,116,10,46,112,115,32,92,110,40,46,83,10,46,102,116,32,92,110,40,68,70,10,0,0,0,0,0,47,112,117,111,114,57,47,53,0,0,0,0,0,0,0,0,47,112,117,111,114,57,47,52,0,0,0,0,0,0,0,0,111,0,0,0,0,0,0,0,84,72,79,82,78,0,0,0,35,102,56,102,56,102,56,0,47,112,117,111,114,57,47,51,0,0,0,0,0,0,0,0,104,114,101,102,0,0,0,0,47,112,117,111,114,57,47,50,0,0,0,0,0,0,0,0,47,112,117,111,114,57,47,49,0,0,0,0,0,0,0,0,47,112,117,111,114,56,47,56,0,0,0,0,0,0,0,0,47,98,114,98,103,49,49,47,49,48,0,0,0,0,0,0,47,112,117,111,114,56,47,55,0,0,0,0,0,0,0,0,47,112,117,111,114,56,47,54,0,0,0,0,0,0,0,0,47,112,117,111,114,56,47,53,0,0,0,0,0,0,0,0,74,117,110,101,0,0,0,0,47,112,117,111,114,56,47,52,0,0,0,0,0,0,0,0,60,33,45,45,32,105,110,115,101,114,116,32,97,110,121,32,111,116,104,101,114,32,104,116,109,108,32,99,111,110,116,101,110,116,32,104,101,114,101,32,45,45,62,10,0,0,0,0,109,97,103,101,110,116,97,0,47,72,101,108,118,101,116,105,99,97,45,66,111,108,100,79,98,108,105,113,117,101,32,115,116,97,114,110,101,116,73,83,79,32,100,101,102,0,0,0,37,115,37,115,32,117,110,115,117,112,112,111,114,116,101,100,10,0,0,0,0,0,0,0,47,112,117,111,114,56,47,51,0,0,0,0,0,0,0,0,47,112,117,111,114,56,47,50,0,0,0,0,0,0,0,0,83,105,103,109,97,0,0,0,35,49,48,49,48,49,48,0,47,112,117,111,114,56,47,49,0,0,0,0,0,0,0,0,108,97,98,101,108,102,108,111,97,116,0,0,0,0,0,0,32,45,45,32,0,0,0,0,47,112,117,111,114,55,47,55,0,0,0,0,0,0,0,0,47,112,117,111,114,55,47,54,0,0,0,0,0,0,0,0,47,112,117,111,114,55,47,53,0,0,0,0,0,0,0,0,47,98,114,98,103,49,49,47,49,0,0,0,0,0,0,0,47,112,117,111,114,55,47,52,0,0,0,0,0,0,0,0,47,112,117,111,114,55,47,51,0,0,0,0,0,0,0,0,47,112,117,111,114,55,47,50,0,0,0,0,0,0,0,0,47,98,114,98,103,49,48,47,57,0,0,0,0,0,0,0,77,97,121,0,0,0,0,0,47,112,117,111,114,55,47,49,0,0,0,0,0,0,0,0,60,68,73,86,32,105,100,61,39,95,86,77,76,50,95,39,32,115,116,121,108,101,61,34,112,111,115,105,116,105,111,110,58,114,101,108,97,116,105,118,101,59,118,105,115,105,98,105,108,105,116,121,58,104,105,100,100,101,110,34,62,10,0,0,108,105,110,101,110,0,0,0,47,72,101,108,118,101,116,105,99,97,45,66,111,108,100,32,115,116,97,114,110,101,116,73,83,79,32,100,101,102,0,0,100,101,102,105,110,101,32,97,116,116,114,115,48,32,37,37,32,37,37,59,32,100,101,102,105,110,101,32,117,110,102,105,108,108,101,100,32,37,37,32,37,37,59,32,100,101,102,105,110,101,32,114,111,117,110,100,101,100,32,37,37,32,37,37,59,32,100,101,102,105,110,101,32,100,105,97,103,111,110,97,108,115,32,37,37,32,37,37,10,0,0,0,0,0,0,0,47,112,117,111,114,54,47,54,0,0,0,0,0,0,0,0,47,97,99,99,101,110,116,52,47,52,0,0,0,0,0,0,47,112,117,111,114,54,47,53,0,0,0,0,0,0,0,0,83,99,97,114,111,110,0,0,35,102,48,102,48,102,48,0,47,112,117,111,114,54,47,52,0,0,0,0,0,0,0,0,32,45,62,32,0,0,0,0,47,112,117,111,114,54,47,51,0,0,0,0,0,0,0,0,118,101,101,0,0,0,0,0,47,112,117,111,114,54,47,50,0,0,0,0,0,0,0,0,47,112,117,111,114,54,47,49,0,0,0,0,0,0,0,0,117,110,115,105,103,110,101,100,32,105,110,116,0,0,0,0,47,112,117,111,114,53,47,53,0,0,0,0,0,0,0,0,47,112,117,111,114,53,47,52,0,0,0,0,0,0,0,0,118,109,108,0,0,0,0,0,110,111,100,101,0,0,0,0,47,112,117,111,114,53,47,51,0,0,0,0,0,0,0,0,47,98,114,98,103,49,48,47,56,0,0,0,0,0,0,0,65,112,114,105,108,0,0,0,47,112,117,111,114,53,47,50,0,0,0,0,0,0,0,0,60,47,68,73,86,62,10,0,108,105,109,101,103,114,101,101,110,0,0,0,0,0,0,0,47,72,101,108,118,101,116,105,99,97,45,79,98,108,105,113,117,101,32,115,116,97,114,110,101,116,73,83,79,32,100,101,102,0,0,0,0,0,0,0,109,105,99,97,0,0,0,0,68,111,116,58,32,91,10,0,47,112,117,111,114,53,47,49,0,0,0,0,0,0,0,0,98,105,115,113,117,101,0,0,47,112,117,111,114,52,47,52,0,0,0,0,0,0,0,0,37,108,102,37,50,115,0,0,82,104,111,0,0,0,0,0,35,101,48,101,48,101,48,0,47,112,117,111,114,52,47,51,0,0,0,0,0,0,0,0,105,110,32,101,100,103,101,32,37,115,37,115,37,115,10,0,32,115,101,116,108,105,110,101,119,105,100,116,104,10,0,0,47,112,117,111,114,52,47,50,0,0,0,0,0,0,0,0,47,112,117,111,114,52,47,49,0,0,0,0,0,0,0,0,47,112,117,111,114,51,47,51,0,0,0,0,0,0,0,0,101,108,108,105,112,115,101,32,97,116,116,114,115,37,100,32,37,115,119,105,100,32,37,46,53,102,32,104,116,32,37,46,53,102,32,97,116,32,40,37,46,53,102,44,37,46,53,102,41,59,10,0,0,0,0,0,99,105,114,99,108,101,32,37,115,32,37,100,44,37,100,44,37,100,10,0,0,0,0,0,47,112,117,111,114,51,47,50,0,0,0,0,0,0,0,0,112,108,97,105,110,58,100,111,116,0,0,0,0,0,0,0,47,112,117,111,114,51,47,49,0,0,0,0,0,0,0,0,106,112,101,58,102,105,103,0,110,111,112,49,0,0,0,0,65,103,114,97,112,104,105,110,102,111,95,116,0,0,0,0,47,112,117,111,114,49,49,47,57,0,0,0,0,0,0,0,106,100,105,97,103,32,62,61,32,48,0,0,0,0,0,0,47,112,117,111,114,49,49,47,56,0,0,0,0,0,0,0,98,97,115,105,99,95,115,116,114,105,110,103,0,0,0,0,60,47,118,58,103,114,111,117,112,62,10,0,0,0,0,0,47,72,101,108,118,101,116,105,99,97,32,115,116,97,114,110,101,116,73,83,79,32,100,101,102,0,0,0,0,0,0,0,109,101,100,95,112,117,114,112,108,101,0,0,0,0,0,0,109,97,120,112,115,104,116,32,61,32,37,102,10,109,97,120,112,115,119,105,100,32,61,32,37,102,10,0,0,0,0,0,47,112,117,111,114,49,49,47,55,0,0,0,0,0,0,0,60,98,117,105,108,116,105,110,62,0,0,0,0,0,0,0,47,112,117,111,114,49,49,47,54,0,0,0,0,0,0,0,115,112,108,105,110,101,115,0,110,0,0,0,0,0,0,0,80,115,105,0,0,0,0,0,35,101,56,101,56,101,56,0,47,112,117,111,114,49,49,47,53,0,0,0,0,0,0,0,118,101,114,116,105,99,101,115,0,0,0,0,0,0,0,0,77,97,114,99,104,0,0,0,47,112,117,111,114,49,49,47,52,0,0,0,0,0,0,0,118,103,0,0,0,0,0,0,47,112,117,111,114,49,49,47,51,0,0,0,0,0,0,0,85,110,104,97,110,100,108,101,100,32,97,100,106,117,115,116,32,111,112,116,105,111,110,32,37,115,10,0,0,0,0,0,75,80,95,68,111,119,110,0,47,112,117,111,114,49,49,47,50,0,0,0,0,0,0,0,47,98,114,98,103,49,48,47,55,0,0,0,0,0,0,0,127,116,111,112,0,0,0,0,110,115,108,105,109,105,116,0,78,68,95,111,114,100,101,114,40,118,41,32,60,32,78,68,95,111,114,100,101,114,40,119,41,0,0,0,0,0,0,0,47,112,117,111,114,49,49,47,49,49,0,0,0,0,0,0,111,114,100,101,114,0,0,0,37,115,32,45,62,32,37,115,58,32,116,97,105,108,32,110,111,116,32,105,110,115,105,100,101,32,116,97,105,108,32,99,108,117,115,116,101,114,32,37,115,10,0,0,0,0,0,0,47,112,117,111,114,49,49,47,49,48,0,0,0,0,0,0,47,112,117,111,114,49,49,47,49,0,0,0,0,0,0,0,99,97,110,110,111,116,32,114,101,97,108,108,111,99,32,100,113,46,112,110,108,115,0,0,70,101,98,114,117,97,114,121,0,0,0,0,0,0,0,0,47,112,117,111,114,49,48,47,57,0,0,0,0,0,0,0,62,10,0,0,0,0,0,0,108,105,103,104,116,121,101,108,108,111,119,0,0,0,0,0,47,84,105,109,101,115,45,66,111,108,100,73,116,97,108,105,99,32,115,116,97,114,110,101,116,73,83,79,32,100,101,102,0,0,0,0,0,0,0,0,109,101,100,105,117,109,119,111,111,100,0,0,0,0,0,0,37,115,32,109,97,120,112,115,104,116,32,97,110,100,32,109,97,120,112,115,119,105,100,32,97,114,101,32,112,114,101,100,101,102,105,110,101,100,32,116,111,32,49,49,46,48,32,97,110,100,32,56,46,53,32,105,110,32,103,112,105,99,10,0,47,112,117,111,114,49,48,47,56,0,0,0,0,0,0,0,47,112,117,111,114,49,48,47,55,0,0,0,0,0,0,0,108,105,98,100,105,114,32,61,32,34,37,115,34,10,0,0,80,114,105,109,101,0,0,0,35,51,48,51,48,51,48,0,47,112,117,111,114,49,48,47,54,0,0,0,0,0,0,0,99,111,109,109,101,110,116,0,98,111,116,104,0,0,0,0,47,112,117,111,114,49,48,47,53,0,0,0,0,0,0,0,85,110,97,98,108,101,32,116,111,32,114,101,99,108,97,105,109,32,98,111,120,32,115,112,97,99,101,32,105,110,32,115,112,108,105,110,101,32,114,111,117,116,105,110,103,32,102,111,114,32,101,100,103,101,32,34,37,115,34,32,45,62,32,34,37,115,34,46,32,83,111,109,101,116,104,105,110,103,32,105,115,32,112,114,111,98,97,98,108,121,32,115,101,114,105,111,117,115,108,121,32,119,114,111,110,103,46,10,0,0,0,0,102,97,105,108,117,114,101,32,109,97,108,108,111,99,39,105,110,103,32,102,111,114,32,114,101,115,117,108,116,32,115,116,114,105,110,103,0,0,0,0,69,78,68,0,0,0,0,0,47,112,117,111,114,49,48,47,52,0,0,0,0,0,0,0,97,100,100,95,116,114,101,101,95,101,100,103,101,58,32,101,109,112,116,121,32,105,110,101,100,103,101,32,108,105,115,116,10,0,0,0,0,0,0,0,47,112,117,111,114,49,48,47,51,0,0,0,0,0,0,0,47,98,114,98,103,49,48,47,54,0,0,0,0,0,0,0,47,112,117,111,114,49,48,47,50,0,0,0,0,0,0,0,44,10,0,0,0,0,0,0,47,112,117,111,114,49,48,47,49,48,0,0,0,0,0,0,47,112,117,111,114,49,48,47,49,0,0,0,0,0,0,0,74,97,110,117,97,114,121,0,47,112,117,98,117,103,110,57,47,57,0,0,0,0,0,0,32,116,97,114,103,101,116,61,34,37,115,34,0,0,0,0,108,105,103,104,116,115,116,101,101,108,98,108,117,101,0,0,47,84,105,109,101,115,45,66,111,108,100,32,115,116,97,114,110,101,116,73,83,79,32,100,101,102,0,0,0,0,0,0,99,108,117,115,116,0,0,0,37,115,32,109,97,120,112,115,104,116,32,97,110,100,32,109,97,120,112,115,119,105,100,32,104,97,118,101,32,110,111,32,109,101,97,110,105,110,103,32,105,110,32,68,87,66,32,50,46,48,44,32,115,101,116,32,112,97,103,101,32,98,111,117,110,100,97,114,105,101,115,32,105,110,32,103,112,105,99,32,97,110,100,32,105,110,32,49,48,116,104,32,69,100,105,116,105,111,110,10,0,0,0,0,47,112,117,98,117,103,110,57,47,56,0,0,0,0,0,0,47,112,117,98,117,103,110,57,47,55,0,0,0,0,0,0,35,102,99,102,99,102,99,0,32,37,115,32,105,110,32,108,105,110,101,32,37,100,32,110,101,97,114,32,39,37,115,39,10,0,0,0,0,0,0,0,47,112,117,98,117,103,110,57,47,54,0,0,0,0,0,0,103,114,111,117,112,0,0,0,98,97,99,107,0,0,0,0,47,112,117,98,117,103,110,57,47,53,0,0,0,0,0,0,47,112,117,98,117,103,110,57,47,52,0,0,0,0,0,0,47,112,117,98,117,103,110,57,47,51,0,0,0,0,0,0,47,98,114,98,103,49,48,47,53,0,0,0,0,0,0,0,47,112,117,98,117,103,110,57,47,50,0,0,0,0,0,0,47,112,117,98,117,103,110,57,47,49,0,0,0,0,0,0,47,112,117,98,117,103,110,56,47,56,0,0,0,0,0,0,47,112,117,98,117,103,110,56,47,55,0,0,0,0,0,0,32,116,105,116,108,101,61,34,37,115,34,0,0,0,0,0,108,105,103,104,116,115,108,97,116,101,103,114,101,121,0,0,47,84,105,109,101,115,45,73,116,97,108,105,99,32,115,116,97,114,110,101,116,73,83,79,32,100,101,102,0,0,0,0,37,115,32,46,80,83,32,119,47,111,32,97,114,103,115,32,99,97,117,115,101,115,32,71,78,85,32,112,105,99,32,116,111,32,115,99,97,108,101,32,100,114,97,119,105,110,103,32,116,111,32,102,105,116,32,56,46,53,120,49,49,32,112,97,112,101,114,59,32,68,87,66,32,100,111,101,115,32,110,111,116,10,0,0,0,0,0,0,47,112,117,98,117,103,110,56,47,54,0,0,0,0,0,0,47,112,117,98,117,103,110,56,47,53,0,0,0,0,0,0,80,104,105,0,0,0,0,0,35,56,48,56,48,56,48,0,47,112,117,98,117,103,110,56,47,52,0,0,0,0,0,0,108,97,121,101,114,0,0,0,102,111,114,119,97,114,100,0,47,112,117,98,117,103,110,56,47,51,0,0,0,0,0,0,47,112,117,98,117,103,110,56,47,50,0,0,0,0,0,0,47,112,117,98,117,103,110,56,47,49,0,0,0,0,0,0,47,98,114,98,103,49,48,47,52,0,0,0,0,0,0,0,47,112,117,98,117,103,110,55,47,55,0,0,0,0,0,0,47,112,117,98,117,103,110,55,47,54,0,0,0,0,0,0,47,112,117,98,117,103,110,55,47,53,0,0,0,0,0,0,68,0,0,0,101,0,0,0,99,0,0,0,0,0,0,0,47,112,117,98,117,103,110,55,47,52,0,0,0,0,0,0,32,104,114,101,102,61,34,37,115,34,0,0,0,0,0,0,108,105,103,104,116,115,108,97,116,101,103,114,97,121,0,0,47,84,105,109,101,115,45,82,111,109,97,110,32,115,116,97,114,110,101,116,73,83,79,32,100,101,102,0,0,0,0,0,108,105,110,101,116,104,105,99,107,32,61,32,48,59,32,111,108,100,108,105,110,101,116,104,105,99,107,32,61,32,108,105,110,101,116,104,105,99,107,10,0,0,0,0,0,0,0,0,47,112,117,98,117,103,110,55,47,51,0,0,0,0,0,0,47,112,117,98,117,103,110,55,47,50,0,0,0,0,0,0,79,117,109,108,0,0,0,0,119,101,100,103,101,100,0,0,47,112,117,98,117,103,110,55,47,49,0,0,0,0,0,0,110,111,106,117,115,116,105,102,121,0,0,0,0,0,0,0,47,112,117,98,117,103,110,54,47,54,0,0,0,0,0,0,47,112,117,98,117,103,110,54,47,53,0,0,0,0,0,0,47,112,117,98,117,103,110,54,47,52,0,0,0,0,0,0,47,98,114,98,103,49,48,47,51,0,0,0,0,0,0,0,103,118,99,86,101,114,115,105,111,110,0,0,0,0,0,0,47,112,117,98,117,103,110,54,47,51,0,0,0,0,0,0,47,112,117,98,117,103,110,54,47,50,0,0,0,0,0,0,47,112,117,98,117,103,110,54,47,49,0,0,0,0,0,0,78,0,0,0,111,0,0,0,118,0,0,0,0,0,0,0,47,112,117,98,117,103,110,53,47,53,0,0,0,0,0,0,60,97,0,0,0,0,0,0,108,105,103,104,116,115,107,121,98,108,117,101,0,0,0,0,125,32,100,101,102,0,0,0,37,115,32,71,78,85,32,112,105,99,32,115,117,112,112,111,114,116,115,32,97,32,108,105,110,101,116,104,105,99,107,32,118,97,114,105,97,98,108,101,32,116,111,32,115,101,116,32,108,105,110,101,32,116,104,105,99,107,110,101,115,115,59,32,68,87,66,32,97,110,100,32,49,48,116,104,32,69,100,46,32,100,111,32,110,111,116,10,0,0,0,0,0,0,0,0,47,112,117,98,117,103,110,53,47,52,0,0,0,0,0,0,47,112,117,98,117,103,110,53,47,51,0,0,0,0,0,0,47,112,117,98,117,103,110,53,47,50,0,0,0,0,0,0,79,116,105,108,100,101,0,0,115,116,114,105,112,101,100,0,105,109,97,103,101,115,99,97,108,101,0,0,0,0,0,0,116,97,112,101,114,101,100,0,47,112,117,98,117,103,110,53,47,49,0,0,0,0,0,0,47,112,117,98,117,103,110,52,47,52,0,0,0,0,0,0,47,112,117,98,117,103,110,52,47,51,0,0,0,0,0,0,47,98,114,98,103,49,48,47,50,0,0,0,0,0,0,0,47,112,117,98,117,103,110,52,47,50,0,0,0,0,0,0,47,112,117,98,117,103,110,52,47,49,0,0,0,0,0,0,47,112,117,98,117,103,110,51,47,51,0,0,0,0,0,0,79,0,0,0,99,0,0,0])
.concat([116,0,0,0,0,0,0,0,47,112,117,98,117,103,110,51,47,50,0,0,0,0,0,0,60,47,97,62,10,0,0,0,108,105,103,104,116,115,101,97,103,114,101,101,110,0,0,0,32,32,32,32,32,32,32,32,99,117,114,114,101,110,116,100,105,99,116,32,101,110,100,32,100,101,102,105,110,101,102,111,110,116,0,0,0,0,0,0,98,111,120,114,97,100,32,61,32,48,32,37,115,32,110,111,32,114,111,117,110,100,101,100,32,99,111,114,110,101,114,115,32,105,110,32,103,114,97,112,104,118,105,122,10,0,0,0,47,112,117,98,117,103,110,51,47,49,0,0,0,0,0,0,47,112,117,98,117,57,47,57,0,0,0,0,0,0,0,0,79,115,108,97,115,104,0,0,114,97,100,105,97,108,0,0,47,112,117,98,117,57,47,56,0,0,0,0,0,0,0,0,102,105,120,101,100,115,105,122,101,0,0,0,0,0,0,0,37,115,45,37,115,0,0,0,47,112,117,98,117,57,47,55,0,0,0,0,0,0,0,0,47,112,117,98,117,57,47,54,0,0,0,0,0,0,0,0,47,112,117,98,117,57,47,53,0,0,0,0,0,0,0,0,47,98,114,98,103,49,48,47,49,48,0,0,0,0,0,0,47,112,117,98,117,57,47,52,0,0,0,0,0,0,0,0,47,112,117,98,117,57,47,51,0,0,0,0,0,0,0,0,47,112,117,98,117,57,47,50,0,0,0,0,0,0,0,0,83,0,0,0,101,0,0,0,112,0,0,0,0,0,0,0,47,112,117,98,117,57,47,49,0,0,0,0,0,0,0,0,60,47,118,58,114,101,99,116,62,10,0,0,0,0,0,0,108,105,103,104,116,115,97,108,109,111,110,0,0,0,0,0,32,32,32,32,32,32,32,32,47,69,110,99,111,100,105,110,103,32,69,110,99,111,100,105,110,103,86,101,99,116,111,114,32,100,101,102,0,0,0,0,37,115,32,71,78,85,32,112,105,99,32,115,117,112,112,111,114,116,115,32,97,32,98,111,120,114,97,100,32,118,97,114,105,97,98,108,101,32,116,111,32,100,114,97,119,32,98,111,120,101,115,32,119,105,116,104,32,114,111,117,110,100,101,100,32,99,111,114,110,101,114,115,59,32,68,87,66,32,97,110,100,32,49,48,116,104,32,69,100,46,32,100,111,32,110,111,116,10,0,0,0,0,0,0,47,112,117,98,117,56,47,56,0,0,0,0,0,0,0,0,34,32,119,105,100,116,104,61,34,37,103,112,120,34,32,104,101,105,103,104,116,61,34,37,103,112,120,34,32,112,114,101,115,101,114,118,101,65,115,112,101,99,116,82,97,116,105,111,61,34,120,77,105,110,89,77,105,110,32,109,101,101,116,34,32,120,61,34,37,103,34,32,121,61,34,37,103,34,0,0,47,112,117,98,117,56,47,55,0,0,0,0,0,0,0,0,79,109,105,99,114,111,110,0,47,112,117,98,117,56,47,54,0,0,0,0,0,0,0,0,100,105,115,116,111,114,116,105,111,110,0,0,0,0,0,0,47,112,117,98,117,56,47,53,0,0,0,0,0,0,0,0,47,112,117,98,117,56,47,52,0,0,0,0,0,0,0,0,47,112,117,98,117,56,47,51,0,0,0,0,0,0,0,0,47,98,114,98,103,49,48,47,49,0,0,0,0,0,0,0,47,112,117,98,117,56,47,50,0,0,0,0,0,0,0,0,47,112,117,98,117,56,47,49,0,0,0,0,0,0,0,0,47,112,117,98,117,55,47,55,0,0,0,0,0,0,0,0,65,0,0,0,117,0,0,0,103,0,0,0,0,0,0,0,47,112,117,98,117,55,47,54,0,0,0,0,0,0,0,0,60,47,99,101,110,116,101,114,62,60,47,118,58,116,101,120,116,98,111,120,62,10,0,0,108,105,103,104,116,112,105,110,107,0,0,0,0,0,0,0,32,32,32,32,32,32,32,32,125,32,102,111,114,97,108,108,0,0,0,0,0,0,0,0,109,101,100,105,117,109,103,111,108,100,101,110,114,111,100,0,97,114,114,111,119,104,101,97,100,32,61,32,55,32,37,115,32,110,111,116,32,117,115,101,100,32,98,121,32,103,114,97,112,104,118,105,122,10,0,0,47,112,117,98,117,55,47,53,0,0,0,0,0,0,0,0,32,116,114,97,110,115,102,111,114,109,61,34,114,111,116,97,116,101,40,37,100,32,37,103,32,37,103,41,34,0,0,0,47,112,117,98,117,55,47,52,0,0,0,0,0,0,0,0,47,97,99,99,101,110,116,52,47,51,0,0,0,0,0,0,79,109,101,103,97,0,0,0,100,105,97,103,111,110,97,108,115,0,0,0,0,0,0,0,47,112,117,98,117,55,47,51,0,0,0,0,0,0,0,0,115,107,101,119,0,0,0,0,47,112,117,98,117,55,47,50,0,0,0,0,0,0,0,0,105,110,118,0,0,0,0,0,47,112,117,98,117,55,47,49,0,0,0,0,0,0,0,0,47,112,117,98,117,54,47,54,0,0,0,0,0,0,0,0,47,98,108,117,101,115,57,47,57,0,0,0,0,0,0,0,105,110,116,0,0,0,0,0,47,112,117,98,117,54,47,53,0,0,0,0,0,0,0,0,47,112,117,98,117,54,47,52,0,0,0,0,0,0,0,0,69,114,114,111,114,32,100,117,114,105,110,103,32,99,111,110,118,101,114,115,105,111,110,32,116,111,32,34,85,84,70,45,56,34,46,32,32,81,117,105,116,105,110,103,46,10,0,0,32,45,100,97,115,104,32,50,0,0,0,0,0,0,0,0,47,112,117,98,117,54,47,51,0,0,0,0,0,0,0,0,74,0,0,0,117,0,0,0,108,0,0,0,0,0,0,0,47,112,117,98,117,54,47,50,0,0,0,0,0,0,0,0,34,62,60,99,101,110,116,101,114,62,0,0,0,0,0,0,108,105,103,104,116,103,114,101,121,0,0,0,0,0,0,0,32,32,32,32,32,32,32,32,123,32,49,32,105,110,100,101,120,32,47,70,73,68,32,110,101,32,123,32,100,101,102,32,125,123,32,112,111,112,32,112,111,112,32,125,32,105,102,101,108,115,101,0,0,0,0,0,109,101,100,105,117,109,102,111,114,101,115,116,103,114,101,101,110,0,0,0,0,0,0,0,37,115,32,97,114,114,111,119,104,101,97,100,32,105,115,32,117,110,100,101,102,105,110,101,100,32,105,110,32,68,87,66,32,50,44,32,105,110,105,116,105,97,108,108,121,32,49,32,105,110,32,103,112,105,99,44,32,50,32,105,110,32,49,48,116,104,32,69,100,105,116,105,111,110,10,0,0,0,0,0,47,112,117,98,117,54,47,49,0,0,0,0,0,0,0,0,34,32,119,105,100,116,104,61,34,37,103,112,120,34,32,104,101,105,103,104,116,61,34,37,103,112,120,34,32,112,114,101,115,101,114,118,101,65,115,112,101,99,116,82,97,116,105,111,61,34,120,77,105,100,89,77,105,100,32,109,101,101,116,34,32,120,61,34,37,103,34,32,121,61,34,37,103,34,0,0,98,101,105,103,101,0,0,0,47,112,117,98,117,53,47,53,0,0,0,0,0,0,0,0,79,103,114,97,118,101,0,0,114,111,117,110,100,101,100,0,47,112,117,98,117,53,47,52,0,0,0,0,0,0,0,0,112,101,114,105,112,104,101,114,105,101,115,0,0,0,0,0,32,93,32,32,37,100,32,102,97,108,115,101,32,37,115,10,0,0,0,0,0,0,0,0,47,112,117,98,117,53,47,51,0,0,0,0,0,0,0,0,47,112,117,98,117,53,47,50,0,0,0,0,0,0,0,0,47,112,117,98,117,53,47,49,0,0,0,0,0,0,0,0,47,98,108,117,101,115,57,47,56,0,0,0,0,0,0,0,32,37,100,0,0,0,0,0,114,101,99,116,32,37,115,32,37,100,44,37,100,32,37,100,44,37,100,10,0,0,0,0,47,112,117,98,117,52,47,52,0,0,0,0,0,0,0,0,99,97,110,111,110,58,100,111,116,0,0,0,0,0,0,0,47,112,117,98,117,52,47,51,0,0,0,0,0,0,0,0,106,112,101,103,58,102,105,103,0,0,0,0,0,0,0,0,110,111,112,0,0,0,0,0,99,108,117,115,116,101,114,0,47,112,117,98,117,52,47,50,0,0,0,0,0,0,0,0,98,101,97,117,116,105,102,121,0,0,0,0,0,0,0,0,74,0,0,0,117,0,0,0,110,0,0,0,0,0,0,0,105,100,105,97,103,32,62,61,32,48,0,0,0,0,0,0,47,112,117,98,117,52,47,49,0,0,0,0,0,0,0,0,99,111,108,111,114,58,35,37,48,50,120,37,48,50,120,37,48,50,120,59,0,0,0,0,108,105,103,104,116,103,114,101,101,110,0,0,0,0,0,0,32,32,32,32,32,32,32,32,100,117,112,32,100,117,112,32,102,105,110,100,102,111,110,116,32,100,117,112,32,108,101,110,103,116,104,32,100,105,99,116,32,98,101,103,105,110,0,0,37,115,32,97,114,114,111,119,104,101,97,100,32,104,97,115,32,110,111,32,109,101,97,110,105,110,103,32,105,110,32,68,87,66,32,50,44,32,97,114,114,111,119,104,101,97,100,32,61,32,55,32,109,97,107,101,115,32,102,105,108,108,101,100,32,97,114,114,111,119,104,101,97,100,115,32,105,110,32,103,112,105,99,32,97,110,100,32,105,110,32,49,48,116,104,32,69,100,105,116,105,111,110,10,0,0,0,0,0,0,0,0,47,112,117,98,117,51,47,51,0,0,0,0,0,0,0,0,60,105,109,97,103,101,32,120,108,105,110,107,58,104,114,101,102,61,34,0,0,0,0,0,65,99,116,105,118,97,116,101,100,32,112,108,117,103,105,110,32,108,105,98,114,97,114,121,58,32,37,115,10,0,0,0,47,112,117,98,117,51,47,50,0,0,0,0,0,0,0,0,79,99,105,114,99,0,0,0,47,112,117,98,117,51,47,49,0,0,0,0,0,0,0,0,115,105,100,101,115,0,0,0,47,112,114,103,110,57,47,57,0,0,0,0,0,0,0,0,100,101,108,116,97,32,60,61,32,48,120,70,70,70,70,0,47,112,114,103,110,57,47,56,0,0,0,0,0,0,0,0,47,112,114,103,110,57,47,55,0,0,0,0,0,0,0,0,47,98,108,117,101,115,57,47,55,0,0,0,0,0,0,0,68,111,119,110,0,0,0,0,127,114,111,111,116,0,0,0,69,100,103,101,32,108,101,110,103,116,104,32,37,102,32,108,97,114,103,101,114,32,116,104,97,110,32,109,97,120,105,109,117,109,32,37,117,32,97,108,108,111,119,101,100,46,10,67,104,101,99,107,32,102,111,114,32,111,118,101,114,119,105,100,101,32,110,111,100,101,40,115,41,46,10,0,0,0,0,0,115,117,114,112,114,105,115,101,10,0,0,0,0,0,0,0,47,112,114,103,110,57,47,54,0,0,0,0,0,0,0,0,69,68,95,116,111,95,118,105,114,116,40,101,41,32,61,61,32,78,85,76,76,0,0,0,49,48,48,48,48,0,0,0,98,101,122,45,62,101,102,108,97,103,0,0,0,0,0,0,47,112,114,103,110,57,47,53,0,0,0,0,0,0,0,0,47,112,114,103,110,57,47,52,0,0,0,0,0,0,0,0,99,97,110,110,111,116,32,109,97,108,108,111,99,32,100,113,46,112,110,108,115,0,0,0,77,0,0,0,97,0,0,0,121,0,0,0,0,0,0,0,47,112,114,103,110,57,47,51,0,0,0,0,0,0,0,0,99,111,108,111,114,58,37,115,59,0,0,0,0,0,0,0,108,105,103,104,116,103,114,97,121,0,0,0,0,0,0,0,47,115,116,97,114,110,101,116,73,83,79,32,123,0,0,0,88,32,101,108,115,101,32,90,10,9,100,101,102,105,110,101,32,115,101,116,102,105,108,108,118,97,108,32,89,32,102,105,108,108,118,97,108,32,61,32,89,59,10,9,100,101,102,105,110,101,32,98,111,108,100,32,89,32,89,59,10,9,100,101,102,105,110,101,32,102,105,108,108,101,100,32,89,32,102,105,108,108,32,89,59,10,90,10,0,0,0,0,0,0,0,0,47,112,114,103,110,57,47,50,0,0,0,0,0,0,0,0,32,37,100,32,37,100,32,37,100,32,37,100,32,37,100,32,37,100,32,37,100,32,37,100,32,37,100,32,37,100,10,0,47,112,114,103,110,57,47,49,0,0,0,0,0,0,0,0,47,103,114,97,112,104,118,105,122,0,0,0,0,0,0,0,79,97,99,117,116,101,0,0,110,111,0,0,0,0,0,0,47,112,114,103,110,56,47,56,0,0,0,0,0,0,0,0,112,101,110,119,105,100,116,104,0,0,0,0,0,0,0,0,47,112,114,103,110,56,47,55,0,0,0,0,0,0,0,0,105,110,32,114,111,117,116,101,115,112,108,105,110,101,115,44,32,80,114,111,117,116,101,115,112,108,105,110,101,32,102,97,105,108,101,100,10,0,0,0,76,97,121,111,117,116,32,119,97,115,32,110,111,116,32,100,111,110,101,10,0,0,0,0,66,69,71,73,78,0,0,0,47,112,114,103,110,56,47,54,0,0,0,0,0,0,0,0,97,100,100,95,116,114,101,101,95,101,100,103,101,58,32,101,109,112,116,121,32,111,117,116,101,100,103,101,32,108,105,115,116,10,0,0,0,0,0,0,47,112,114,103,110,56,47,53,0,0,0,0,0,0,0,0,47,98,108,117,101,115,57,47,54,0,0,0,0,0,0,0,47,112,114,103,110,56,47,52,0,0,0,0,0,0,0,0,32,91,0,0,0,0,0,0,47,112,114,103,110,56,47,51,0,0,0,0,0,0,0,0,47,112,114,103,110,56,47,50,0,0,0,0,0,0,0,0,65,0,0,0,112,0,0,0,114,0,0,0,0,0,0,0,47,112,114,103,110,56,47,49,0,0,0,0,0,0,0,0,32,102,111,110,116,45,115,105,122,101,58,32,37,46,50,102,112,116,59,0,0,0,0,0,108,105,103,104,116,103,111,108,100,101,110,114,111,100,121,101,108,108,111,119,0,0,0,0,37,32,83,101,116,32,117,112,32,73,83,79,32,76,97,116,105,110,32,49,32,99,104,97,114,97,99,116,101,114,32,101,110,99,111,100,105,110,103,0,9,37,115,9,115,111,114,114,121,44,32,116,104,101,32,103,114,111,102,102,32,102,111,108,107,115,32,99,104,97,110,103,101,100,32,103,112,105,99,59,32,115,101,110,100,32,97,110,121,32,99,111,109,112,108,97,105,110,116,32,116,111,32,116,104,101,109,59,10,0,0,0,47,112,114,103,110,55,47,55,0,0,0,0,0,0,0,0,37,100,32,37,100,32,37,100,32,37,100,32,37,100,32,37,100,32,37,100,32,37,100,32,37,100,32,37,46,49,102,32,37,100,32,37,100,32,37,100,32,37,100,32,37,100,32,37,100,10,32,37,100,32,37,115,10,0,0,0,0,0,0,0,47,112,114,103,110,55,47,54,0,0,0,0,0,0,0,0,79,69,108,105,103,0,0,0,37,115,58,37,100,58,32,37,115,32,105,110,32,108,105,110,101,32,37,100,32,110,101,97,114,32,39,37,115,39,10,0,47,112,114,103,110,55,47,53,0,0,0,0,0,0,0,0,120,108,97,98,101,108,0,0,47,112,114,103,110,55,47,52,0,0,0,0,0,0,0,0,47,112,114,103,110,55,47,51,0,0,0,0,0,0,0,0,47,112,114,103,110,55,47,50,0,0,0,0,0,0,0,0,47,98,108,117,101,115,57,47,53,0,0,0,0,0,0,0,47,112,114,103,110,55,47,49,0,0,0,0,0,0,0,0,47,112,114,103,110,54,47,54,0,0,0,0,0,0,0,0,47,112,114,103,110,54,47,53,0,0,0,0,0,0,0,0,77,0,0,0,97,0,0,0,114,0,0,0,0,0,0,0,47,112,114,103,110,54,47,52,0,0,0,0,0,0,0,0,102,111,110,116,45,115,116,121,108,101,58,32,37,115,59,0,108,105,103,104,116,99,121,97,110,0,0,0,0,0,0,0,69,110,99,111,100,105,110,103,86,101,99,116,111,114,32,52,53,32,47,104,121,112,104,101,110,32,112,117,116,0,0,0,109,97,110,100,97,114,105,110,111,114,97,110,103,101,0,0,9,37,115,9,105,110,115,116,97,108,108,32,97,32,109,111,114,101,32,114,101,99,101,110,116,32,118,101,114,115,105,111,110,32,111,102,32,103,112,105,99,32,111,114,32,115,119,105,116,99,104,32,116,111,32,68,87,66,32,111,114,32,49,48,116,104,32,69,100,105,116,105,111,110,32,112,105,99,59,10,0,0,0,0,0,0,0,0,47,112,114,103,110,54,47,51,0,0,0,0,0,0,0,0,125,10,0,0,0,0,0,0,47,112,114,103,110,54,47,50,0,0,0,0,0,0,0,0,78,117,0,0,0,0,0,0,110,111,100,101,32,37,115,44,32,112,111,114,116,32,37,115,32,117,110,114,101,99,111,103,110,105,122,101,100,10,0,0,47,112,114,103,110,54,47,49,0,0,0,0,0,0,0,0,102,111,110,116,99,111,108,111,114,0,0,0,0,0,0,0,47,112,114,103,110,53,47,53,0,0,0,0,0,0,0,0,47,112,114,103,110,53,47,52,0,0,0,0,0,0,0,0,47,112,114,103,110,53,47,51,0,0,0,0,0,0,0,0,47,98,108,117,101,115,57,47,52,0,0,0,0,0,0,0,47,112,114,103,110,53,47,50,0,0,0,0,0,0,0,0,47,112,114,103,110,53,47,49,0,0,0,0,0,0,0,0,47,112,114,103,110,52,47,52,0,0,0,0,0,0,0,0,70,0,0,0,101,0,0,0,98,0,0,0,0,0,0,0,47,112,114,103,110,52,47,51,0,0,0,0,0,0,0,0,102,111,110,116,45,115,116,114,101,116,99,104,58,32,37,115,59,0,0,0,0,0,0,0,108,105,103,104,116,99,111,114,97,108,0,0,0,0,0,0,73,83,79,76,97,116,105,110,49,69,110,99,111,100,105,110,103,32,48,32,50,53,53,32,103,101,116,105,110,116,101,114,118,97,108,32,112,117,116,105,110,116,101,114,118,97,108,0,9,37,115,32,105,102,32,121,111,117,32,117,115,101,32,103,112,105,99,32,97,110,100,32,105,116,32,98,97,114,102,115,32,111,110,32,101,110,99,111,117,110,116,101,114,105,110,103,32,34,115,111,108,105,100,34,44,10,0,0,0,0,0,0,47,112,114,103,110,52,47,50,0,0,0,0,0,0,0,0,32,32,125,10,0,0,0,0,47,112,114,103,110,52,47,49,0,0,0,0,0,0,0,0,78,116,105,108,100,101,0,0,110,111,100,101,32,37,115,44,32,112,111,114,116,32,37,115,44,32,117,110,114,101,99,111,103,110,105,122,101,100,32,99,111,109,112,97,115,115,32,112,111,105,110,116,32,39,37,115,39,32,45,32,105,103,110,111,114,101,100,10,0,0,0,0,47,112,114,103,110,51,47,51,0,0,0,0,0,0,0,0,102,111,110,116,110,97,109,101,0,0,0,0,0,0,0,0,47,112,114,103,110,51,47,50,0,0,0,0,0,0,0,0,47,112,114,103,110,51,47,49,0,0,0,0,0,0,0,0,47,112,114,103,110,49,49,47,57,0,0,0,0,0,0,0,47,98,108,117,101,115,57,47,51,0,0,0,0,0,0,0,65,103,100,101,115,99,95,116,0,0,0,0,0,0,0,0,47,112,114,103,110,49,49,47,56,0,0,0,0,0,0,0,47,112,114,103,110,49,49,47,55,0,0,0,0,0,0,0,47,112,114,103,110,49,49,47,54,0,0,0,0,0,0,0,74,0,0,0,97,0,0,0,110,0,0,0,0,0,0,0,47,112,114,103,110,49,49,47,53,0,0,0,0,0,0,0,102,111,110,116,45,119,101,105,103,104,116,58,32,37,115,59,0,0,0,0,0,0,0,0,108,105,103,104,116,98,108,117,101,0,0,0,0,0,0,0,32,69,110,99,111,100,105,110,103,86,101,99,116,111,114,32,48,0,0,0,0,0,0,0,105,102,32,102,105,108,108,118,97,108,32,62,32,48,46,52,32,116,104,101,110,32,88,10,9,100,101,102,105,110,101,32,115,101,116,102,105,108,108,118,97,108,32,89,32,102,105,108,108,118,97,108,32,61,32,49,32,45,32,89,59,10,9,100,101,102,105,110,101,32,98,111,108,100,32,89,32,116,104,105,99,107,110,101,115,115,32,50,32,89,59,10,0,0,0,0,47,112,114,103,110,49,49,47,52,0,0,0,0,0,0,0,32,32,32,32,116,101,120,116,117,114,101,32,73,109,97,103,101,84,101,120,116,117,114,101,32,123,32,117,114,108,32,34,37,115,34,32,125,10,0,0,47,112,114,103,110,49,49,47,51,0,0,0,0,0,0,0,77,117,0,0,0,0,0,0,95,0,0,0,0,0,0,0,47,112,114,103,110,49,49,47,50,0,0,0,0,0,0,0,102,111,110,116,115,105,122,101,0,0,0,0,0,0,0,0,47,112,114,103,110,49,49,47,49,49,0,0,0,0,0,0,98,122,46,115,105,122,101,32,37,32,51,32,61,61,32,49,0,0,0,0,0,0,0,0,108,105,98,103,114,97,112,104,58,32,112,97,114,115,101,114,32,108,111,115,116,32,114,111,111,116,32,103,114,97,112,104,10,0,0,0,0,0,0,0,47,112,114,103,110,49,49,47,49,48,0,0,0,0,0,0,47,112,114,103,110,49,49,47,49,0,0,0,0,0,0,0,47,98,108,117,101,115,57,47,50,0,0,0,0,0,0,0,47,112,114,103,110,49,48,47,57,0,0,0,0,0,0,0,47,112,114,103,110,49,48,47,56,0,0,0,0,0,0,0,47,112,114,103,110,49,48,47,55,0,0,0,0,0,0,0,47,112,114,103,110,49,48,47,54,0,0,0,0,0,0,0,102,111,110,116,45,102,97,109,105,108,121,58,32,39,37,115,39,59,0,0,0,0,0,0,108,101,109,111,110,99,104,105,102,102,111,110,0,0,0,0,47,69,110,99,111,100,105,110,103,86,101,99,116,111,114,32,50,53,54,32,97,114,114,97,121,32,100,101,102,0,0,0,108,105,103,104,116,119,111,111,100,0,0,0,0,0,0,0,37,115,32,71,78,85,32,112,105,99,32,118,115,46,32,49,48,116,104,32,69,100,105,116,105,111,110,32,100,92,40,101,39,116,101,110,116,101,10,0,47,112,114,103,110,49,48,47,53,0,0,0,0,0,0,0,32,32,32,32,125,10,0,0,47,112,114,103,110,49,48,47,52,0,0,0,0,0,0,0,76,97,109,98,100,97,0,0,37,46,53,103,32,37,46,53,103,32,116,114,97,110,115,108,97,116,101,32,110,101,119,112,97,116,104,32,117,115,101,114,95,115,104,97,112,101,95,37,100,10,0,0,0,0,0,0,47,112,114,103,110,49,48,47,51,0,0,0,0,0,0,0,68,0,0,0,101,0,0,0,99,0,0,0,101,0,0,0,109,0,0,0,98,0,0,0,101,0,0,0,114,0,0,0,0,0,0,0,0,0,0,0,98,122,46,115,105,122,101,32,62,32,48,0,0,0,0,0,47,112,114,103,110,49,48,47,50,0,0,0,0,0,0,0,47,112,114,103,110,49,48,47,49,48,0,0,0,0,0,0,47,112,114,103,110,49,48,47,49,0,0,0,0,0,0,0,47,98,108,117,101,115,57,47,49,0,0,0,0,0,0,0,47,112,105,121,103,57,47,57,0,0,0,0,0,0,0,0,47,112,105,121,103,57,47,56,0,0,0,0,0,0,0,0,47,112,105,121,103,57,47,55,0,0,0,0,0,0,0,0,78,0,0,0,111,0,0,0,118,0,0,0,101,0,0,0,109,0,0,0,98,0,0,0,101,0,0,0,114,0,0,0,0,0,0,0,0,0,0,0,47,112,105,121,103,57,47,54,0,0,0,0,0,0,0,0,60,118,58,116,101,120,116,98,111,120,32,105,110,115,101,116,61,34,48,44,48,44,48,44,48,34,32,115,116,121,108,101,61,34,112,111,115,105,116,105,111,110,58,97,98,115,111,108,117,116,101,59,32,118,45,116,101,120,116,45,119,114,97,112,112,105,110,103,58,39,102,97,108,115,101,39,59,112,97,100,100,105,110,103,58,39,48,39,59,0,0,0,0,0,0,0,108,97,119,110,103,114,101,101,110,0,0,0,0,0,0,0,109,97,114,107,0,0,0,0,114,101,115,101,116,32,37,115,32,115,101,116,32,116,111,32,107,110,111,119,110,32,115,116,97,116,101,10,0,0,0,0,47,112,105,121,103,57,47,53,0,0,0,0,0,0,0,0,32,32,32,32,32,32,32,32,100,105,102,102,117,115,101,67,111,108,111,114,32,49,32,49,32,49,10,0,0,0,0,0,47,112,105,121,103,57,47,52,0,0,0,0,0,0,0,0,75,97,112,112,97,0,0,0,77,114,101,99,111,114,100,0,47,112,105,121,103,57,47,51,0,0,0,0,0,0,0,0,102,105,108,108,99,111,108,111,114,0,0,0,0,0,0,0,115,112,108,45,62,115,105,122,101,32,62,32,48,0,0,0,47,112,105,121,103,57,47,50,0,0,0,0,0,0,0,0,47,112,105,121,103,57,47,49,0,0,0,0,0,0,0,0,47,112,105,121,103,56,47,56,0,0,0,0,0,0,0,0,47,98,108,117,101,115,56,47,56,0,0,0,0,0,0,0,47,112,105,121,103,56,47,55,0,0,0,0,0,0,0,0,47,112,105,121,103,56,47,54,0,0,0,0,0,0,0,0,47,112,105,121,103,56,47,53,0,0,0,0,0,0,0,0,79,0,0,0,99,0,0,0,116,0,0,0,111,0,0,0,98,0,0,0,101,0,0,0,114,0,0,0,0,0,0,0,47,112,105,121,103,56,47,52,0,0,0,0,0,0,0,0,32,115,116,114,111,107,101,100,61,34,102,97,108,115,101,34,32,102,105,108,108,101,100,61,34,102,97,108,115,101,34,62,10,0,0,0,0,0,0,0,108,97,118,101,110,100,101,114,98,108,117,115,104,0,0,0,47,115,101,116,117,112,76,97,116,105,110,49,32,123,0,0,108,105,103,104,116,95,112,117,114,112,108,101,0,0,0,0,105,102,32,98,111,120,114,97,100,32,62,32,49,46,48,32,38,38,32,100,97,115,104,119,105,100,32,60,32,48,46,48,55,53,32,116,104,101,110,32,88,10,9,102,105,108,108,118,97,108,32,61,32,49,59,10,9,100,101,102,105,110,101,32,102,105,108,108,32,89,32,89,59,10,9,100,101,102,105,110,101,32,115,111,108,105,100,32,89,32,89,59,10,9,100,101,102,105,110,101,32,114,101,115,101,116,32,89,32,115,99,97,108,101,61,49,46,48,32,89,59,10,88,10,0,0,0,0,47,112,105,121,103,56,47,51,0,0,0,0,0,0,0,0,32,32,32,32,32,32,97,109,98,105,101,110,116,73,110,116,101,110,115,105,116,121,32,48,46,51,51,10,0,0,0,0,37,108,102,44,37,108,102,44,37,108,102,44,37,108,102,0,47,112,105,121,103,56,47,50,0,0,0,0,0,0,0,0,47,97,99,99,101,110,116,52,47,50,0,0,0,0,0,0,73,117,109,108,0,0,0,0,114,101,99,111,114,100,0,0,47,112,105,121,103,56,47,49,0,0,0,0,0,0,0,0,99,111,108,111,114,0,0,0,115,101,116,108,105,110,101,119,105,100,116,104,0,49,0,0,47,112,105,121,103,55,47,55,0,0,0,0,0,0,0,0,47,112,105,121,103,55,47,54,0,0,0,0,0,0,0,0,47,112,105,121,103,55,47,53,0,0,0,0,0,0,0,0,47,98,108,117,101,115,56,47,55,0,0,0,0,0,0,0,117,110,115,105,103,110,101,100,32,115,104,111,114,116,0,0,47,112,105,121,103,55,47,52,0,0,0,0,0,0,0,0,47,112,105,121,103,55,47,51,0,0,0,0,0,0,0,0,38,35,51,57,59,0,0,0,32,45,100,97,115,104,32,53,0,0,0,0,0,0,0,0,47,112,105,121,103,55,47,50,0,0,0,0,0,0,0,0,47,112,105,121,103,55,47,49,0,0,0,0,0,0,0,0,60,118,58,114,101,99,116,32,115,116,121,108,101,61,34,112,111,115,105,116,105,111,110,58,97,98,115,111,108,117,116,101,59,32,0,0,0,0,0,0,108,97,118,101,110,100,101,114,0,0,0,0,0,0,0,0,37,115,32,68,87,66,32,50,32,99,111,109,112,97,116,105,98,105,108,105,116,121,32,100,101,102,105,110,105,116,105,111,110,115,10,0,0,0,0,0,47,112,105,121,103,54,47,54,0,0,0,0,0,0,0,0,32,32,32,32,109,97,116,101,114,105,97,108,32,77,97,116,101,114,105,97,108,32,123,10,0,0,0,0,0,0,0,0,97,122,117,114,101,0,0,0,47,112,105,121,103,54,47,53,0,0,0,0,0,0,0,0,99,97,110,110,111,116,32,99,111,109,112,105,108,101,32,114,101,103,117,108,97,114,32,101,120,112,114,101,115,115,105,111,110,32,37,115,0,0,0,0,73,111,116,97,0,0,0,0,108,112,114,111,109,111,116,101,114,0,0,0,0,0,0,0,47,112,105,121,103,54,47,52,0,0,0,0,0,0,0,0,83,0,0,0,101,0,0,0,112,0,0,0,116,0,0,0,101,0,0,0,109,0,0,0,98,0,0,0,101,0,0,0,114,0,0,0,0,0,0,0,115,111,108,105,100,0,0,0,32,93,32,32,37,100,32,116,114,117,101,32,37,115,10,0,47,112,105,121,103,54,47,51,0,0,0,0,0,0,0,0,47,112,105,121,103,54,47,50,0,0,0,0,0,0,0,0,47,112,105,121,103,54,47,49,0,0,0,0,0,0,0,0,47,98,108,117,101,115,56,47,54,0,0,0,0,0,0,0,32,37,115,10,0,0,0,0,99,109,97,112,120,95,110,112,58,109,97,112,0,0,0,0,47,112,105,121,103,53,47,53,0,0,0,0,0,0,0,0,103,118,58,100,111,116,0,0,47,112,105,121,103,53,47,52,0,0,0,0,0,0,0,0,103,105,102,58,102,105,103,0,111,115,97,103,101,0,0,0,71,114,97,112,104,32,37,115,32,104,97,115,32,97,114,114,97,121,32,112,97,99,107,105,110,103,32,119,105,116,104,32,117,115,101,114,32,118,97,108,117,101,115,32,98,117,116,32,110,111,32,34,115,111,114,116,118,34,32,97,116,116,114,105,98,117,116,101,115,32,97,114,101,32,100,101,102,105,110,101,100,46,0,0,0,0,0,0,47,112,105,121,103,53,47,51,0,0,0,0,0,0,0,0,113,117,97,100,116,114,101,101,0,0,0,0,0,0,0,0,65,0,0,0,117,0,0,0,103,0,0,0,117,0,0,0,115,0,0,0,116,0,0,0,0,0,0,0,0,0,0,0,47,112,105,121,103,53,47,50,0,0,0,0,0,0,0,0,83,112,97,114,115,101,77,97,116,114,105,120,95,105,115,95,115,121,109,109,101,116,114,105,99,40,65,44,32,70,65,76,83,69,41,32,38,38,32,65,45,62,116,121,112,101,32,61,61,32,77,65,84,82,73,88,95,84,89,80,69,95,82,69,65,76,0,0,0,0,0,0,60,47,118,58,111,118,97,108,62,10,0,0,0,0,0,0,107,104,97,107,105,0,0,0,68,111,116,68,105,99,116,32,98,101,103,105,110,0,0,0,37,115,32,114,101,115,101,116,32,119,111,114,107,115,32,105,110,32,103,112,105,99,32,97,110,100,32,49,48,116,104,32,101,100,105,116,105,111,110,44,32,98,117,116,32,105,115,110,39,116,32,100,101,102,105,110,101,100,32,105,110,32,68,87,66,32,50,10,0,0,0,0,115,116,97,114,116,61,37,115,32,110,111,116,32,115,117,112,112,111,114,116,101,100,32,119,105,116,104,32,109,111,100,101,61,115,101,108,102,32,45,32,105,103,110,111,114,101,100,10,0,0,0,0,0,0,0,0,47,112,105,121,103,53,47,49,0,0,0,0,0,0,0,0,32,32,97,112,112,101,97,114,97,110,99,101,32,65,112,112,101,97,114,97,110,99,101,32,123,10,0,0,0,0,0,0,102,97,105,108,101,100,32,116,111,32,114,101,115,111,108,118,101,32,37,115,32,105,110,32,37,115,10,0,0,0,0,0,116,97,105,108,95,108,112,0,47,112,105,121,103,52,47,52,0,0,0,0,0,0,0,0,73,103,114,97,118,101,0,0,114,97,114,114,111,119,0,0,98,98,0,0,0,0,0,0,47,112,105,121,103,52,47,51,0,0,0,0,0,0,0,0,47,112,105,121,103,52,47,50,0,0,0,0,0,0,0,0,47,112,105,121,103,52,47,49,0,0,0,0,0,0,0,0,87,97,114,110,105,110,103,58,32,110,111,100,101,32,37,115,44,32,112,111,115,105,116,105,111,110,32,37,115,44,32,101,120,112,101,99,116,101,100,32,116,119,111,32,102,108,111,97,116,115,10,0,0,0,0,0,47,112,105,121,103,51,47,51,0,0,0,0,0,0,0,0,47,98,108,117,101,115,56,47,53,0,0,0,0,0,0,0,75,80,95,85,112,0,0,0,115,101,97,114,99,104,115,105,122,101,0,0,0,0,0,0,99,111,110,116,97,105,110,95,110,111,100,101,115,32,99,108,117,115,116,32,37,115,32,114,97,110,107,32,37,100,32,109,105,115,115,105,110,103,32,110,111,100,101,10,0,0,0,0,105,110,115,116,97,108,108,95,105,110,95,114,97,110,107,44,32,108,105,110,101,32,37,100,58,32,71,68,95,114,97,110,107,40,103,41,91,37,100,93,46,118,32,43,32,78,68,95,111,114,100,101,114,40,37,115,41,32,91,37,100,93,32,62,32,71,68,95,114,97,110,107,40,103,41,91,37,100,93,46,97,118,32,43,32,71,68,95,114,97,110,107,40,82,111,111,116,41,91,37,100,93,46,97,110,32,91,37,100,93,10,0,47,112,105,121,103,51,47,50,0,0,0,0,0,0,0,0,109,101,114,103,101,95,111,110,101,119,97,121,32,103,108,105,116,99,104,10,0,0,0,0,78,111,32,108,105,98,122,32,115,117,112,112,111,114,116,10,0,0,0,0,0,0,0,0,98,101,122,45,62,115,102,108,97,103,0,0,0,0,0,0,47,112,105,121,103,51,47,49,0,0,0,0,0,0,0,0,47,112,105,121,103,49,49,47,57,0,0,0,0,0,0,0,99,97,110,110,111,116,32,114,101,97,108,108,111,99,32,111,112,115,0,0,0,0,0,0,74,0,0,0,117,0,0,0,108,0,0,0,121,0,0,0,0,0,0,0,0,0,0,0,47,112,105,121,103,49,49,47,56,0,0,0,0,0,0,0,32,119,105,100,116,104,58,32,37,46,50,102,59,32,104,101,105,103,104,116,58,32,37,46,50,102,34,0,0,0,0,0,105,118,111,114,121,0,0,0,47,68,111,116,68,105,99,116,32,50,48,48,32,100,105,99,116,32,100,101,102,0,0,0,37,115,32,68,87,66,32,50,32,100,111,101,115,110,39,116,32,117,115,101,32,102,105,108,108,32,97,110,100,32,100,111,101,115,110,39,116,32,100,101,102,105,110,101,32,102,105,108,108,118,97,108,10,0,0,0,47,112,105,121,103,49,49,47,55,0,0,0,0,0,0,0,83,104,97,112,101,32,123,10,0,0,0,0,0,0,0,0,104,101,97,100,95,108,112,0,47,112,105,121,103,49,49,47,54,0,0,0,0,0,0,0,47,46,108,105,98,115,0,0,73,99,105,114,99,0,0,0,108,97,114,114,111,119,0,0,47,112,105,121,103,49,49,47,53,0,0,0,0,0,0,0,47,112,105,121,103,49,49,47,52,0,0,0,0,0,0,0,35,37,50,120,37,50,120,37,50,120,37,50,120,0,0,0,105,110,32,114,111,117,116,101,115,112,108,105,110,101,115,44,32,80,115,104,111,114,116,101,115,116,112,97,116,104,32,102,97,105,108,101,100,10,0,0,69,79,70,0,0,0,0,0,47,112,105,121,103,49,49,47,51,0,0,0,0,0,0,0,97,100,100,95,116,114,101,101,95,101,100,103,101,58,32,109,105,115,115,105,110,103,32,116,114,101,101,32,101,100,103,101,10,0,0,0,0,0,0,0,107,105,110,100,32,61,61,32,76,84,95,78,79,78,69,0,47,112,105,121,103,49,49,47,50,0,0,0,0,0,0,0,47,98,108,117,101,115,56,47,52,0,0,0,0,0,0,0,103,114,97,112,104,0,0,0,47,112,105,121,103,49,49,47,49,49,0,0,0,0,0,0,93,0,0,0,0,0,0,0,47,112,105,121,103,49,49,47,49,48,0,0,0,0,0,0,47,112,105,121,103,49,49,47,49,0,0,0,0,0,0,0,74,0,0,0,117,0,0,0,110,0,0,0,101,0,0,0,0,0,0,0,0,0,0,0,47,112,105,121,103,49,48,47,57,0,0,0,0,0,0,0,32,108,101,102,116,58,32,37,46,50,102,59,32,116,111,112,58,32,37,46,50,102,59,0,105,110,100,105,103,111,0,0,37,37,66,101,103,105,110,80,114,111,108,111,103,0,0,0,104,117,110,116,101,114,115,103,114,101,101,110,0,0,0,0,37,115,32,102,105,108,108,32,104,97,115,32,110,111,32,109,101,97,110,105,110,103,32,105,110,32,68,87,66,32,50,44,32,103,112,105,99,32,99,97,110,32,117,115,101,32,102,105,108,108,32,111,114,32,102,105,108,108,101,100,44,32,49,48,116,104,32,69,100,105,116,105,111,110,32,117,115,101,115,32,102,105,108,108,32,111,110,108,121,10,0,0,0,0,0,0,47,112,105,121,103,49,48,47,56,0,0,0,0,0,0,0,115,121,110,116,97,120,32,101,114,114,111,114,32,105,110,32,112,111,115,32,97,116,116,114,105,98,117,116,101,32,102,111,114,32,101,100,103,101,32,40,37,115,44,37,115,41,10,0,47,112,105,121,103,49,48,47,55,0,0,0,0,0,0,0,73,97,99,117,116,101,0,0,114,112,114,111,109,111,116,101,114,0,0,0,0,0,0,0,47,112,105,121,103,49,48,47,54,0,0,0,0,0,0,0,68,105,110,103,98,97,116,115,0,0,0,0,0,0,0,0,47,112,105,121,103,49,48,47,53,0,0,0,0,0,0,0,47,112,105,121,103,49,48,47,52,0,0,0,0,0,0,0,47,112,105,121,103,49,48,47,51,0,0,0,0,0,0,0,47,98,108,117,101,115,56,47,51,0,0,0,0,0,0,0,47,112,105,121,103,49,48,47,50,0,0,0,0,0,0,0,47,112,105,121,103,49,48,47,49,48,0,0,0,0,0,0,47,112,105,121,103,49,48,47,49,0,0,0,0,0,0,0,47,98,108,117,101,115,56,47,50,0,0,0,0,0,0,0,47,112,97,115,116,101,108,50,56,47,56,0,0,0,0,0,32,32,60,118,58,111,118,97,108,32,115,116,121,108,101,61,34,112,111,115,105,116,105,111,110,58,97,98,115,111,108,117,116,101,59,0,0,0,0,0,105,110,100,105,97,110,114,101,100,0,0,0,0,0,0,0,91,32,123,67,97,116,97,108,111,103,125,32,60,60,32,47,85,82,73,32,60,60,32,47,66,97,115,101,32,40,37,115,41,32,62,62,32,62,62,10,47,80,85,84,32,112,100,102,109,97,114,107,10,0,0,0,115,116,100,58,58,98,97,100,95,99,97,115,116,0,0,0,37,115,32,102,105,108,108,118,97,108,32,105,115,32,48,46,51,32,105,110,32,49,48,116,104,32,69,100,105,116,105,111,110,32,40,102,105,108,108,32,48,32,109,101,97,110,115,32,98,108,97,99,107,41,44,32,48,46,53,32,105,110,32,103,112,105,99,32,40,102,105,108,108,32,48,32,109,101,97,110,115,32,119,104,105,116,101,41,44,32,117,110,100,101,102,105,110,101,100,32,105,110,32,68,87,66,32,50,10,0,0,0,47,112,97,115,116,101,108,50,56,47,55,0,0,0,0,0,111,98,106,0,0,0,0,0,37,108,102,44,37,108,102,37,110,0,0,0,0,0,0,0,47,112,97,115,116,101,108,50,56,47,54,0,0,0,0,0,71,97,109,109,97,0,0,0,115,105,103,110,97,116,117,114,101,0,0,0,0,0,0,0,47,112,97,115,116,101,108,50,56,47,53,0,0,0,0,0,103,114,97,100,105,101,110,116,97,110,103,108,101,0,0,0,112,97,103,101,100,105,114,0,47,112,97,115,116,101,108,50,56,47,52,0,0,0,0,0,47,112,97,115,116,101,108,50,56,47,51,0,0,0,0,0,47,112,97,115,116,101,108,50,56,47,50,0,0,0,0,0,97,103,100,105,99,116,111,102,58,32,117,110,107,110,111,119,110,32,107,105,110,100,32,37,100,10,0,0,0,0,0,0,47,112,97,115,116,101,108,50,56,47,49,0,0,0,0,0,47,112,97,115,116,101,108,50,55,47,55,0,0,0,0,0,47,112,97,115,116,101,108,50,55,47,54,0,0,0,0,0,47,98,108,117,101,115,56,47,49,0,0,0,0,0,0,0,65,0,0,0,112,0,0,0,114,0,0,0,105,0,0,0,108,0,0,0,0,0,0,0,47,112,97,115,116,101,108,50,55,47,53,0,0,0,0,0,120,32,101,32,34,47,62,0,104,111,116,112,105,110,107,0,115,101,116,117,112,76,97,116,105,110,49,10,0,0,0,0,103,114,101,101,110,99,111,112,112,101,114,0,0,0,0,0,37,115,32,100,97,115,104,119,105,100,32,105,115,32,48,46,49,32,105,110,32,49,48,116,104,32,69,100,105,116,105,111,110,44,32,48,46,48,53,32,105,110,32,68,87,66,32,50,32,97,110,100,32,105,110,32,103,112,105,99,10,0,0,0,47,112,97,115,116,101,108,50,55,47,52,0,0,0,0,0,112,111,115,32,97,116,116,114,105,98,117,116,101,32,102,111,114,32,101,100,103,101,32,40,37,115,44,37,115,41,32,100,111,101,115,110,39,116,32,104,97,118,101,32,51,110,43,49,32,112,111,105,110,116,115,10,0,0,0,0,0,0,0,0,47,112,97,115,116,101,108,50,55,47,51,0,0,0,0,0,69,117,109,108,0,0,0,0,97,115,115,101,109,98,108,121,0,0,0,0,0,0,0,0,47,112,97,115,116,101,108,50,55,47,50,0,0,0,0,0,111,114,100,101,114,105,110,103,0,0,0,0,0,0,0,0,109,101,100,105,117,109,0,0,66,76,0,0,0,0,0,0,47,112,97,115,116,101,108,50,55,47,49,0,0,0,0,0,47,112,97,115,116,101,108,50,54,47,54,0,0,0,0,0,47,112,97,115,116,101,108,50,54,47,53,0,0,0,0,0,97,103,97,112,112,108,121,58,32,117,110,107,110,111,119,110,32,111,98,106,101,99,116,32,116,121,112,101,32,37,100,10,0,0,0,0,0,0,0,0,65,103,101,100,103,101,95,116,0,0,0,0,0,0,0,0,47,112,97,115,116,101,108,50,54,47,52,0,0,0,0,0,47,112,97,115,116,101,108,50,54,47,51,0,0,0,0,0,47,112,97,115,116,101,108,50,54,47,50,0,0,0,0,0,77,0,0,0,97,0,0,0,114,0,0,0,99,0,0,0,104,0,0,0,0,0,0,0,47,112,97,115,116,101,108,50,54,47,49,0,0,0,0,0,108,32,0,0,0,0,0,0,104,111,110,101,121,100,101,119,0,0,0,0,0,0,0,0,37,37,69,110,100,67,111,109,109,101,110,116,115,10,115,97,118,101,10,0,0,0,0,0,37,115,32,98,111,120,114,97,100,32,105,115,32,110,111,119,32,48,46,48,32,105,110,32,103,112,105,99,44,32,101,108,115,101,32,105,116,32,114,101,109,97,105,110,115,32,50,46,48,10,0,0,0,0,0,0,47,112,97,115,116,101,108,50,53,47,53,0,0,0,0,0,117,115,101,114,95,115,104,97,112,101,95,37,100,10,0,0,32,101,44,37,108,102,44,37,108,102,37,110,0,0,0,0,47,112,97,115,116,101,108,50,53,47,52,0,0,0,0,0,69,116,97,0,0,0,0,0,110,111,118,101,114,104,97,110,103,0,0,0,0,0,0,0,47,112,97,115,116,101,108,50,53,47,51,0,0,0,0,0,114,101,115,111,108,117,116,105,111,110,0,0,0,0,0,0,47,112,97,115,116,101,108,50,53,47,50,0,0,0,0,0,85,82,87,32,67,104,97,110,99,101,114,121,32,76,0,0,112,97,100,0,0,0,0,0,47,112,97,115,116,101,108,50,53,47,49,0,0,0,0,0,118,111,105,100,0,0,0,0,47,112,97,115,116,101,108,50,52,47,52,0,0,0,0,0,47,98,108,117,101,115,55,47,55,0,0,0,0,0,0,0,47,112,97,115,116,101,108,50,52,47,51,0,0,0,0,0,47,112,97,115,116,101,108,50,52,47,50,0,0,0,0,0,47,112,97,115,116,101,108,50,52,47,49,0,0,0,0,0,70,0,0,0,101,0,0,0,98,0,0,0,114,0,0,0])
.concat([117,0,0,0,97,0,0,0,114,0,0,0,121,0,0,0,0,0,0,0,0,0,0,0,47,112,97,115,116,101,108,50,51,47,51,0,0,0,0,0,37,46,48,102,32,37,46,48,102,32,0,0,0,0,0,0,103,114,101,121,0,0,0,0,37,37,37,37,66,111,117,110,100,105,110,103,66,111,120,58,32,37,100,32,37,100,32,37,100,32,37,100,10,0,0,0,103,114,97,121,57,53,0,0,115,99,97,108,101,61,49,46,48,32,37,115,32,114,101,113,117,105,114,101,100,32,102,111,114,32,99,111,109,112,97,114,105,115,111,110,115,10,0,0,47,112,97,115,116,101,108,50,51,47,50,0,0,0,0,0,103,115,97,118,101,32,37,103,32,37,103,32,116,114,97,110,115,108,97,116,101,32,110,101,119,112,97,116,104,10,0,0,115,44,37,108,102,44,37,108,102,37,110,0,0,0,0,0,47,112,97,115,116,101,108,50,51,47,49,0,0,0,0,0,69,112,115,105,108,111,110,0,116,104,114,101,101,112,111,118,101,114,104,97,110,103,0,0,47,112,97,115,116,101,108,49,57,47,57,0,0,0,0,0,100,112,105,0,0,0,0,0,90,97,112,102,67,104,97,110,99,101,114,121,45,77,101,100,105,117,109,73,116,97,108,105,99,0,0,0,0,0,0,0,47,112,97,115,116,101,108,49,57,47,56,0,0,0,0,0,47,112,97,115,116,101,108,49,57,47,55,0,0,0,0,0,47,112,97,115,116,101,108,49,57,47,54,0,0,0,0,0,47,98,108,117,101,115,55,47,54,0,0,0,0,0,0,0,47,112,97,115,116,101,108,49,57,47,53,0,0,0,0,0,47,112,97,115,116,101,108,49,57,47,52,0,0,0,0,0,47,112,97,115,116,101,108,49,57,47,51,0,0,0,0,0,74,0,0,0,97,0,0,0,110,0,0,0,117,0,0,0,97,0,0,0,114,0,0,0,121,0,0,0,0,0,0,0,47,112,97,115,116,101,108,49,57,47,50,0,0,0,0,0,38,97,109,112,59,0,0,0,32,102,105,108,108,101,100,61,34,102,97,108,115,101,34,32,0,0,0,0,0,0,0,0,103,114,101,101,110,121,101,108,108,111,119,0,0,0,0,0,37,37,66,111,117,110,100,105,110,103,66,111,120,58,32,40,97,116,101,110,100,41,10,0,103,114,97,121,57,48,0,0,98,111,120,114,97,100,61,50,46,48,32,37,115,32,119,105,108,108,32,98,101,32,114,101,115,101,116,32,116,111,32,48,46,48,32,98,121,32,103,112,105,99,32,111,110,108,121,10,0,0,0,0,0,0,0,0,47,112,97,115,116,101,108,49,57,47,49,0,0,0,0,0,93,32,32,37,100,32,102,97,108,115,101,32,37,115,10,0,47,112,97,115,116,101,108,49,56,47,56,0,0,0,0,0,69,103,114,97,118,101,0,0,102,105,118,101,112,111,118,101,114,104,97,110,103,0,0,0,47,112,97,115,116,101,108,49,56,47,55,0,0,0,0,0,99,111,110,99,101,110,116,114,97,116,101,0,0,0,0,0,47,112,97,115,116,101,108,49,56,47,54,0,0,0,0,0,116,107,0,0,0,0,0,0,47,112,97,115,116,101,108,49,56,47,53,0,0,0,0,0,47,112,97,115,116,101,108,49,56,47,52,0,0,0,0,0,47,98,108,117,101,115,55,47,53,0,0,0,0,0,0,0,47,112,97,115,116,101,108,49,56,47,51,0,0,0,0,0,47,112,97,115,116,101,108,49,56,47,50,0,0,0,0,0,47,112,97,115,116,101,108,49,56,47,49,0,0,0,0,0,47,112,97,115,116,101,108,49,55,47,55,0,0,0,0,0,34,32,0,0,0,0,0,0,37,37,80,97,103,101,115,58,32,49,10,0,0,0,0,0,103,114,97,121,56,53,0,0,37,115,32,110,111,110,45,102,97,116,97,108,32,114,117,110,45,116,105,109,101,32,112,105,99,32,118,101,114,115,105,111,110,32,100,101,116,101,114,109,105,110,97,116,105,111,110,44,32,118,101,114,115,105,111,110,32,50,10,0,0,0,0,0,47,112,97,115,116,101,108,49,55,47,54,0,0,0,0,0,93,32,32,37,100,32,116,114,117,101,32,37,115,10,0,0,112,105,110,0,0,0,0,0,47,112,97,115,116,101,108,49,55,47,53,0,0,0,0,0,47,97,99,99,101,110,116,52,47,49,0,0,0,0,0,0,69,99,105,114,99,0,0,0,114,101,115,116,114,105,99,116,105,111,110,115,105,116,101,0,47,112,97,115,116,101,108,49,55,47,52,0,0,0,0,0,99,108,117,115,116,101,114,114,97,110,107,0,0,0,0,0,84,104,101,32,99,104,97,114,97,99,116,101,114,32,39,37,99,39,32,97,112,112,101,97,114,115,32,105,110,32,98,111,116,104,32,116,104,101,32,108,97,121,101,114,115,101,112,32,97,110,100,32,108,97,121,101,114,108,105,115,116,115,101,112,32,97,116,116,114,105,98,117,116,101,115,32,45,32,108,97,121,101,114,108,105,115,116,115,101,112,32,105,103,110,111,114,101,100,46,10,0,0,0,0,47,112,97,115,116,101,108,49,55,47,51,0,0,0,0,0,47,112,97,115,116,101,108,49,55,47,50,0,0,0,0,0,47,112,97,115,116,101,108,49,55,47,49,0,0,0,0,0,47,98,108,117,101,115,55,47,52,0,0,0,0,0,0,0,115,104,111,114,116,0,0,0,47,112,97,115,116,101,108,49,54,47,54,0,0,0,0,0,47,112,97,115,116,101,108,49,54,47,53,0,0,0,0,0,38,113,117,111,116,59,0,0,32,45,102,105,108,108,32,0,47,112,97,115,116,101,108,49,54,47,52,0,0,0,0,0,47,98,108,117,101,115,55,47,51,0,0,0,0,0,0,0,47,112,97,115,116,101,108,49,54,47,51,0,0,0,0,0,32,102,105,108,108,101,100,61,34,116,114,117,101,34,32,102,105,108,108,99,111,108,111,114,61,34,0,0,0,0,0,0,37,37,80,97,103,101,115,58,32,40,97,116,101,110,100,41,10,0,0,0,0,0,0,0,103,114,97,121,56,48,0,0,37,115,32,100,111,110,39,116,32,99,104,97,110,103,101,32,97,110,121,116,104,105,110,103,32,98,101,108,111,119,32,116,104,105,115,32,108,105,110,101,32,105,110,32,116,104,105,115,32,100,114,97,119,105,110,103,10,0,0,0,0,0,0,0,47,112,97,115,116,101,108,49,54,47,50,0,0,0,0,0,37,103,32,37,103,32,0,0,97,113,117,97,109,97,114,105,110,101,0,0,0,0,0,0,47,112,97,115,116,101,108,49,54,47,49,0,0,0,0,0,34,37,115,34,32,119,97,115,32,110,111,116,32,102,111,117,110,100,32,97,115,32,97,32,102,105,108,101,32,111,114,32,97,115,32,97,32,115,104,97,112,101,32,108,105,98,114,97,114,121,32,109,101,109,98,101,114,10,0,0,0,0,0,0,69,97,99,117,116,101,0,0,112,114,105,109,101,114,115,105,116,101,0,0,0,0,0,0,47,112,97,115,116,101,108,49,53,47,53,0,0,0,0,0,108,97,110,100,115,99,97,112,101,0,0,0,0,0,0,0,80,77,0,0,0,0,0,0,44,0,0,0,0,0,0,0,32,0,0,0,0,0,0,0,47,112,97,115,116,101,108,49,53,47,52,0,0,0,0,0,105,110,118,105,115,105,98,108,101,0,0,0,0,0,0,0,47,112,97,115,116,101,108,49,53,47,51,0,0,0,0,0,98,97,107,101,114,115,99,104,111,99,0,0,0,0,0,0,47,112,97,115,116,101,108,49,53,47,50,0,0,0,0,0,105,109,97,112,95,110,112,58,109,97,112,0,0,0,0,0,47,112,97,115,116,101,108,49,53,47,49,0,0,0,0,0,100,111,116,58,100,111,116,0,47,112,97,115,116,101,108,49,52,47,52,0,0,0,0,0,112,110,103,58,102,105,103,0,112,97,116,99,104,119,111,114,107,0,0,0,0,0,0,0,115,111,114,116,118,0,0,0,47,112,97,115,116,101,108,49,52,47,51,0,0,0,0,0,47,98,108,117,101,115,55,47,50,0,0,0,0,0,0,0,115,118,103,0,0,0,0,0,115,109,111,111,116,104,105,110,103,0,0,0,0,0,0,0,65,77,0,0,0,0,0,0,47,112,97,115,116,101,108,49,52,47,50,0,0,0,0,0,34,0,0,0,0,0,0,0,103,111,108,100,101,110,114,111,100,0,0,0,0,0,0,0,37,37,37,37,84,105,116,108,101,58,32,37,115,10,0,0,103,114,97,121,55,53,0,0,46,110,114,32,83,70,32,37,46,48,102,10,115,99,97,108,101,116,104,105,99,107,110,101,115,115,32,61,32,37,46,48,102,10,0,0,0,0,0,0,47,112,97,115,116,101,108,49,52,47,49,0,0,0,0,0,105,110,118,97,108,105,100,32,112,108,117,103,105,110,32,112,97,116,104,32,34,37,115,34,10,0,0,0,0,0,0,0,47,112,97,115,116,101,108,49,51,47,51,0,0,0,0,0,69,84,72,0,0,0,0,0,112,114,111,116,101,105,110,115,116,97,98,0,0,0,0,0,108,112,0,0,0,0,0,0,47,112,97,115,116,101,108,49,51,47,50,0,0,0,0,0,105,112,0,0,0,0,0,0,111,114,105,101,110,116,97,116,105,111,110,0,0,0,0,0,84,105,109,101,115,0,0,0,108,97,121,101,114,108,105,115,116,115,101,112,0,0,0,0,47,112,97,115,116,101,108,49,51,47,49,0,0,0,0,0,71,114,97,112,104,118,105,122,32,98,117,105,108,116,32,119,105,116,104,111,117,116,32,97,110,121,32,116,114,105,97,110,103,117,108,97,116,105,111,110,32,108,105,98,114,97,114,121,10,0,0,0,0,0,0,0,47,112,97,105,114,101,100,57,47,57,0,0,0,0,0,0,111,118,101,114,108,97,112,0,47,112,97,105,114,101,100,57,47,56,0,0,0,0,0,0,112,116,114,0,0,0,0,0,85,112,0,0,0,0,0,0,108,101,118,101,108,32,103,114,97,112,104,32,114,101,99,0,105,110,115,116,97,108,108,95,105,110,95,114,97,110,107,44,32,108,105,110,101,32,37,100,58,32,114,97,110,107,32,37,100,32,110,111,116,32,105,110,32,114,97,110,107,32,114,97,110,103,101,32,91,37,100,44,37,100,93,10,0,0,0,0,47,112,97,105,114,101,100,57,47,55,0,0,0,0,0,0,102,105,110,100,95,102,97,115,116,95,110,111,100,101,40,103,44,32,110,41,0,0,0,0,118,32,61,61,32,110,0,0,46,46,47,46,46,47,46,46,47,108,105,98,47,100,111,116,103,101,110,47,99,111,109,112,111,117,110,100,46,99,0,0,47,112,97,105,114,101,100,57,47,54,0,0,0,0,0,0,47,112,97,105,114,101,100,57,47,53,0,0,0,0,0,0,99,97,110,110,111,116,32,109,97,108,108,111,99,32,111,112,115,0,0,0,0,0,0,0,47,98,108,117,101,115,55,47,49,0,0,0,0,0,0,0,47,112,97,105,114,101,100,57,47,52,0,0,0,0,0,0,103,111,108,100,0,0,0,0,37,100,32,37,100,32,115,101,116,108,97,121,101,114,10,0,103,114,97,121,55,48,0,0,37,115,32,116,111,32,99,104,97,110,103,101,32,100,114,97,119,105,110,103,32,115,105,122,101,44,32,109,117,108,116,105,112,108,121,32,116,104,101,32,119,105,100,116,104,32,97,110,100,32,104,101,105,103,104,116,32,111,110,32,116,104,101,32,46,80,83,32,108,105,110,101,32,97,98,111,118,101,32,97,110,100,32,116,104,101,32,110,117,109,98,101,114,32,111,110,32,116,104,101,32,116,119,111,32,108,105,110,101,115,32,98,101,108,111,119,32,40,114,111,117,110,100,101,100,32,116,111,32,116,104,101,32,110,101,97,114,101,115,116,32,105,110,116,101,103,101,114,41,32,98,121,32,97,32,115,99,97,108,101,32,102,97,99,116,111,114,10,0,0,0,0,0,0,0,0,47,112,97,105,114,101,100,57,47,51,0,0,0,0,0,0,117,115,45,62,110,97,109,101,0,0,0,0,0,0,0,0,47,112,97,105,114,101,100,57,47,50,0,0,0,0,0,0,68,101,108,116,97,0,0,0,47,108,105,98,103,118,99,46,0,0,0,0,0,0,0,0,112,114,111,116,101,97,115,101,115,105,116,101,0,0,0,0,47,0,0,0,0,0,0,0,47,112,97,105,114,101,100,57,47,49,0,0,0,0,0,0,114,111,116,97,116,101,0,0,58,9,32,0,0,0,0,0,47,112,97,105,114,101,100,56,47,56,0,0,0,0,0,0,99,117,115,116,111,109,0,0,105,110,32,114,111,117,116,101,115,112,108,105,110,101,115,44,32,101,100,103,101,32,105,115,32,97,32,108,111,111,112,32,97,116,32,37,115,10,0,0,99,97,110,39,116,32,111,112,101,110,32,108,105,98,114,97,114,121,32,102,105,108,101,32,37,115,10,0,0,0,0,0,37,100,32,111,117,116,32,111,102,32,37,100,32,101,120,116,101,114,105,111,114,32,108,97,98,101,108,115,32,112,111,115,105,116,105,111,110,101,100,46,10,0,0,0,0,0,0,0,47,112,97,105,114,101,100,56,47,55,0,0,0,0,0,0,47,112,97,105,114,101,100,56,47,54,0,0,0,0,0,0,117,112,100,97,116,101,58,32,109,105,115,109,97,116,99,104,101,100,32,108,99,97,32,105,110,32,116,114,101,101,117,112,100,97,116,101,115,10,0,0,46,46,47,46,46,47,46,46,47,108,105,98,47,99,111,109,109,111,110,47,108,97,98,101,108,115,46,99,0,0,0,0,46,46,47,46,46,47,46,46,47,108,105,98,47,103,118,99,47,103,118,117,115,101,114,115,104,97,112,101,46,99,0,0,99,58,47,101,109,115,99,114,105,112,116,101,110,47,115,121,115,116,101,109,47,105,110,99,108,117,100,101,92,101,109,115,99,114,105,112,116,101,110,47,98,105,110,100,46,104,0,0,47,112,97,105,114,101,100,56,47,53,0,0,0,0,0,0,32,91,107,101,121,61,0,0,47,112,97,105,114,101,100,56,47,52,0,0,0,0,0,0,109,101,109,111,114,121,32,101,120,104,97,117,115,116,101,100,0,0,0,0,0,0,0,0,47,112,97,105,114,101,100,56,47,51,0,0,0,0,0,0,116,101,120,116,108,97,121,111,117,116,0,0,0,0,0,0,47,98,108,117,101,115,54,47,54,0,0,0,0,0,0,0,104,101,108,118,101,116,105,99,97,0,0,0,0,0,0,0,80,0,0,0,77,0,0,0,0,0,0,0,0,0,0,0,47,112,97,105,114,101,100,56,47,50,0,0,0,0,0,0,99,32,0,0,0,0,0,0,103,104,111,115,116,119,104,105,116,101,0,0,0,0,0,0,91,32,47,67,114,111,112,66,111,120,32,91,37,100,32,37,100,32,37,100,32,37,100,93,32,47,80,65,71,69,83,32,112,100,102,109,97,114,107,10,0,0,0,0,0,0,0,0,103,114,97,121,54,53,0,0,46,80,83,32,37,46,53,102,32,37,46,53,102,10,0,0,47,112,97,105,114,101,100,56,47,49,0,0,0,0,0,0,117,115,0,0,0,0,0,0,112,115,0,0,0,0,0,0,73,108,108,101,103,97,108,32,118,97,108,117,101,32,37,115,32,102,111,114,32,97,116,116,114,105,98,117,116,101,32,34,109,111,100,101,34,32,105,110,32,103,114,97,112,104,32,37,115,32,45,32,105,103,110,111,114,101,100,10,0,0,0,0,47,112,97,105,114,101,100,55,47,55,0,0,0,0,0,0,68,97,103,103,101,114,0,0,114,110,97,115,116,97,98,0,47,112,97,105,114,101,100,55,47,54,0,0,0,0,0,0,99,101,110,116,101,114,0,0,102,97,110,116,97,115,121,0,108,97,121,101,114,115,101,112,0,0,0,0,0,0,0,0,47,112,97,105,114,101,100,55,47,53,0,0,0,0,0,0,47,112,97,105,114,101,100,55,47,52,0,0,0,0,0,0,47,112,97,105,114,101,100,55,47,51,0,0,0,0,0,0,97,103,101,100,103,101,115,101,116,0,0,0,0,0,0,0,47,112,97,105,114,101,100,55,47,50,0,0,0,0,0,0,47,112,97,105,114,101,100,55,47,49,0,0,0,0,0,0,47,112,97,105,114,101,100,54,47,54,0,0,0,0,0,0,47,98,108,117,101,115,54,47,53,0,0,0,0,0,0,0,65,0,0,0,77,0,0,0,0,0,0,0,0,0,0,0,47,112,97,105,114,101,100,54,47,53,0,0,0,0,0,0,37,115,37,46,48,102,44,37,46,48,102,32,0,0,0,0,103,97,105,110,115,98,111,114,111,0,0,0,0,0,0,0,99,97,110,118,97,115,32,115,105,122,101,32,40,37,100,44,37,100,41,32,101,120,99,101,101,100,115,32,80,68,70,32,108,105,109,105,116,32,40,37,100,41,10,9,40,115,117,103,103,101,115,116,32,115,101,116,116,105,110,103,32,97,32,98,111,117,110,100,105,110,103,32,98,111,120,32,115,105,122,101,44,32,115,101,101,32,100,111,116,40,49,41,41,10,0,0,103,114,97,121,54,48,0,0,114,111,116,97,116,105,111,110,0,0,0,0,0,0,0,0,47,112,97,105,114,101,100,54,47,52,0,0,0,0,0,0,108,111,97,100,105,109,97,103,101,0,0,0,0,0,0,0,106,111,98,0,0,0,0,0,104,105,101,114,0,0,0,0,47,112,97,105,114,101,100,54,47,51,0,0,0,0,0,0,67,104,105,0,0,0,0,0,114,105,98,111,115,105,116,101,0,0,0,0,0,0,0,0,47,112,97,105,114,101,100,54,47,50,0,0,0,0,0,0,112,97,103,101,0,0,0,0,97,108,108,0,0,0,0,0,47,112,97,105,114,101,100,54,47,49,0,0,0,0,0,0,47,112,97,105,114,101,100,53,47,53,0,0,0,0,0,0,47,112,97,105,114,101,100,53,47,52,0,0,0,0,0,0,97,103,101,100,103,101,103,101,116,0,0,0,0,0,0,0,47,112,97,105,114,101,100,53,47,51,0,0,0,0,0,0,47,112,97,105,114,101,100,53,47,50,0,0,0,0,0,0,47,112,97,105,114,101,100,53,47,49,0,0,0,0,0,0,47,98,108,117,101,115,54,47,52,0,0,0,0,0,0,0,47,112,97,105,114,101,100,52,47,52,0,0,0,0,0,0,109,32,0,0,0,0,0,0,47,112,97,105,114,101,100,52,47,51,0,0,0,0,0,0,37,103,32,37,103,32,115,101,116,95,115,99,97,108,101,32,37,100,32,114,111,116,97,116,101,32,37,103,32,37,103,32,116,114,97,110,115,108,97,116,101,10,0,0,0,0,0,0,103,114,97,121,53,53,0,0,93,10,46,80,69,10,0,0,100,101,118,105,99,101,0,0,46,46,47,46,46,47,46,46,47,112,108,117,103,105,110,47,99,111,114,101,47,103,118,108,111,97,100,105,109,97,103,101,95,99,111,114,101,46,99,0,109,97,106,111,114,0,0,0,47,112,97,105,114,101,100,52,47,50,0,0,0,0,0,0,67,99,101,100,105,108,0,0,105,110,115,117,108,97,116,111,114,0,0,0,0,0,0,0,47,112,97,105,114,101,100,52,47,49,0,0,0,0,0,0,115,105,122,101,0,0,0,0,84,104,101,32,108,97,121,101,114,115,101,108,101,99,116,32,97,116,116,114,105,98,117,116,101,32,34,37,115,34,32,100,111,101,115,32,110,111,116,32,109,97,116,99,104,32,97,110,121,32,108,97,121,101,114,32,115,112,101,99,105,102,101,100,32,98,121,32,116,104,101,32,108,97,121,101,114,115,32,97,116,116,114,105,98,117,116,101,32,45,32,105,103,110,111,114,101,100,46,10,0,0,0,0,47,112,97,105,114,101,100,51,47,51,0,0,0,0,0,0,47,112,97,105,114,101,100,51,47,50,0,0,0,0,0,0,47,112,97,105,114,101,100,51,47,49,0,0,0,0,0,0,97,103,110,111,100,101,115,101,116,0,0,0,0,0,0,0,65,103,110,111,100,101,95,116,0,0,0,0,0,0,0,0,47,112,97,105,114,101,100,49,50,47,57,0,0,0,0,0,47,112,97,105,114,101,100,49,50,47,56,0,0,0,0,0,47,112,97,105,114,101,100,49,50,47,55,0,0,0,0,0,47,98,108,117,101,115,54,47,51,0,0,0,0,0,0,0,112,105,99,0,0,0,0,0,47,112,97,105,114,101,100,49,50,47,54,0,0,0,0,0,47,62,60,47,118,58,115,104,97,112,101,62,10,0,0,0,102,111,114,101,115,116,103,114,101,101,110,0,0,0,0,0,103,115,97,118,101,10,37,100,32,37,100,32,37,100,32,37,100,32,98,111,120,112,114,105,109,32,99,108,105,112,32,110,101,119,112,97,116,104,10,0,103,114,97,121,53,48,0,0,90,97,112,102,68,105,110,103,98,97,116,115,0,0,0,0,47,112,97,105,114,101,100,49,50,47,53,0,0,0,0,0,32,47,62,10,0,0,0,0,75,75,0,0,0,0,0,0,47,112,97,105,114,101,100,49,50,47,52,0,0,0,0,0,66,101,116,97,0,0,0,0,117,116,114,0,0,0,0,0,47,112,97,105,114,101,100,49,50,47,51,0,0,0,0,0,102,111,110,116,110,97,109,101,115,0,0,0,0,0,0,0,108,97,121,101,114,115,101,108,101,99,116,0,0,0,0,0,109,97,112,0,0,0,0,0,47,112,97,105,114,101,100,49,50,47,50,0,0,0,0,0,47,112,97,105,114,101,100,49,50,47,49,50,0,0,0,0,47,112,97,105,114,101,100,49,50,47,49,49,0,0,0,0,97,103,110,111,100,101,103,101,116,0,0,0,0,0,0,0,47,112,97,105,114,101,100,49,50,47,49,48,0,0,0,0,47,112,97,105,114,101,100,49,50,47,49,0,0,0,0,0,102,105,103,0,0,0,0,0,47,112,97,105,114,101,100,49,49,47,57,0,0,0,0,0,47,98,108,117,101,115,54,47,50,0,0,0,0,0,0,0,47,112,97,105,114,101,100,49,49,47,56,0,0,0,0,0,60,118,58,112,97,116,104,32,32,118,61,34,0,0,0,0,102,108,111,114,97,108,119,104,105,116,101,0,0,0,0,0,37,100,32,37,100,32,37,100,32,98,101,103,105,110,112,97,103,101,10,0,0,0,0,0,103,114,97,121,52,53,0,0,83,121,109,98,111,108,0,0,47,112,97,105,114,101,100,49,49,47,55,0,0,0,0,0,108,97,121,111,117,116,0,0,60,118,58,105,109,97,103,101,32,115,114,99,61,34,37,115,34,32,115,116,121,108,101,61,34,32,112,111,115,105,116,105,111,110,58,97,98,115,111,108,117,116,101,59,32,119,105,100,116,104,58,37,46,50,102,59,32,104,101,105,103,104,116,58,37,46,50,102,59,32,108,101,102,116,58,37,46,50,102,32,59,32,116,111,112,58,37,46,50,102,34,0,0,0,0,0,109,111,100,101,0,0,0,0,47,112,97,105,114,101,100,49,49,47,54,0,0,0,0,0,99,111,114,101,0,0,0,0,65,117,109,108,0,0,0,0,116,101,114,109,105,110,97,116,111,114,0,0,0,0,0,0,47,112,97,105,114,101,100,49,49,47,53,0,0,0,0,0,115,104,111,119,98,111,120,101,115,0,0,0,0,0,0,0,108,97,121,101,114,115,0,0,47,112,97,105,114,101,100,49,49,47,52,0,0,0,0,0,47,112,97,105,114,101,100,49,49,47,51,0,0,0,0,0,47,112,97,105,114,101,100,49,49,47,50,0,0,0,0,0,97,103,115,101,116,0,0,0,47,112,97,105,114,101,100,49,49,47,49,49,0,0,0,0,47,112,97,105,114,101,100,49,49,47,49,48,0,0,0,0,47,112,97,105,114,101,100,49,49,47,49,0,0,0,0,0,47,98,108,117,101,115,54,47,49,0,0,0,0,0,0,0,47,112,97,105,114,101,100,49,48,47,57,0,0,0,0,0,32,62,0,0,0,0,0,0,112,110,103,58,115,118,103,0,102,105,114,101,98,114,105,99,107,0,0,0,0,0,0,0,60,60,32,47,80,97,103,101,83,105,122,101,32,91,37,100,32,37,100,93,32,62,62,32,115,101,116,112,97,103,101,100,101,118,105,99,101,10,0,0,110,101,97,116,111,95,108,97,121,111,117,116,0,0,0,0,103,114,97,121,52,48,0,0,84,105,109,101,115,45,82,111,109,97,110,0,0,0,0,0,47,112,97,105,114,101,100,49,48,47,56,0,0,0,0,0,114,101,110,100,101,114,0,0,106,112,103,58,118,109,108,0,85,110,107,110,111,119,110,32,118,97,108,117,101,32,37,115,32,102,111,114,32,97,116,116,114,105,98,117,116,101,32,34,109,111,100,101,108,34,32,105,110,32,103,114,97,112,104,32,37,115,32,45,32,105,103,110,111,114,101,100,10,0,0,0,110,101,97,116,111,0,0,0,112,114,105,115,109,0,0,0,47,112,97,105,114,101,100,49,48,47,55,0,0,0,0,0,100,111,116,95,108,97,121,111,117,116,0,0,0,0,0,0,65,116,105,108,100,101,0,0,99,100,115,0,0,0,0,0,47,112,97,105,114,101,100,49,48,47,54,0,0,0,0,0,101,113,117,97,108,108,121,0,80,97,108,97,116,105,110,111,32,76,105,110,111,116,121,112,101,0,0,0,0,0,0,0,100,103,101,115,102,105,114,115,116,0,0,0,0,0,0,0,47,112,97,105,114,101,100,49,48,47,53,0,0,0,0,0,47,112,97,105,114,101,100,49,48,47,52,0,0,0,0,0,47,112,97,105,114,101,100,49,48,47,51,0,0,0,0,0,97,103,103,101,116,0,0,0,97,114,101,97,0,0,0,0,114,111,111,116,0,0,0,0,47,112,97,105,114,101,100,49,48,47,50,0,0,0,0,0,46,46,47,46,46,47,46,46,47,108,105,98,47,115,102,100,112,103,101,110,47,117,110,105,102,111,114,109,95,115,116,114,101,115,115,46,99,0,0,0,47,112,97,105,114,101,100,49,48,47,49,48,0,0,0,0,47,112,97,105,114,101,100,49,48,47,49,0,0,0,0,0,47,98,108,117,101,115,53,47,53,0,0,0,0,0,0,0,47,111,114,114,100,57,47,57,0,0,0,0,0,0,0,0,32,119,105,100,116,104,58,32,37,100,59,32,104,101,105,103,104,116,58,32,37,100,34,0,100,111,100,103,101,114,98,108,117,101,0,0,0,0,0,0,80,111,114,116,114,97,105,116,0,0,0,0,0,0,0,0,103,114,97,121,51,53,0,0,80,97,108,97,116,105,110,111,45,66,111,108,100,73,116,97,108,105,99,0,0,0,0,0,47,111,114,114,100,57,47,56,0,0,0,0,0,0,0,0,106,112,101,58,118,109,108,0,105,115,32,105,110,97,112,112,114,111,112,114,105,97,116,101,46,32,82,101,118,101,114,116,105,110,103,32,116,111,32,116,104,101,32,115,104,111,114,116,101,115,116,32,112,97,116,104,32,109,111,100,101,108,46,10,0,0,0,0,0,0,0,0,121,120,32,112,115,101,117,100,111,45,111,114,116,104,111,103,111,110,97,108,32,99,111,110,115,116,114,97,105,110,116,115,0,0,0,0,0,0,0,0,47,111,114,114,100,57,47,55,0,0,0,0,0,0,0,0,47,97,99,99,101,110,116,51,47,51,0,0,0,0,0,0,65,114,105,110,103,0,0,0,112,114,111,109,111,116,101,114,0,0,0,0,0,0,0,0,47,111,114,114,100,57,47,54,0,0,0,0,0,0,0,0,95,37,100,0,0,0,0,0,111,100,101,115,102,105,114,115,116,0,0,0,0,0,0,0,47,111,114,114,100,57,47,53,0,0,0,0,0,0,0,0,47,111,114,114,100,57,47,52,0,0,0,0,0,0,0,0,112,111,115,0,0,0,0,0,47,111,114,114,100,57,47,51,0,0,0,0,0,0,0,0,97,103,100,101,108,101,100,103,101,0,0,0,0,0,0,0,47,111,114,114,100,57,47,50,0,0,0,0,0,0,0,0,46,46,47,46,46,47,46,46,47,108,105,98,47,115,102,100,112,103,101,110,47,112,111,115,116,95,112,114,111,99,101,115,115,46,99,0,0,0,0,0,117,110,115,105,103,110,101,100,32,99,104,97,114,0,0,0,47,111,114,114,100,57,47,49,0,0,0,0,0,0,0,0,38,35,49,54,48,59,0,0,32,99,114,101,97,116,101,32,108,105,110,101,32,0,0,0,47,111,114,114,100,56,47,56,0,0,0,0,0,0,0,0,47,98,108,117,101,115,53,47,52,0,0,0,0,0,0,0,46,46,47,46,46,47,46,46,47,108,105,98,47,115,102,100,112,103,101,110,47,81,117,97,100,84,114,101,101,46,99,0,47,111,114,114,100,56,47,55,0,0,0,0,0,0,0,0,48,0,0,0,0,0,0,0,114,105,102,102,0,0,0,0,100,105,109,103,114,101,121,0,76,97,110,100,115,99,97,112,101,0,0,0,0,0,0,0,103,114,97,121,51,48,0,0,80,97,108,97,116,105,110,111,45,73,116,97,108,105,99,0,47,111,114,114,100,56,47,54,0,0,0,0,0,0,0,0,106,112,101,103,58,118,109,108,0,0,0,0,0,0,0,0,46,46,47,46,46,47,46,46,47,108,105,98,47,115,102,100,112,103,101,110,47,80,114,105,111,114,105,116,121,81,117,101,117,101,46,99,0,0,0,0,101,100,103,101,115,32,105,110,32,103,114,97,112,104,32,37,115,32,104,97,118,101,32,110,111,32,108,101,110,32,97,116,116,114,105,98,117,116,101,46,32,72,101,110,99,101,44,32,116,104,101,32,109,100,115,32,109,111,100,101,108,10,0,0,112,111,114,116,104,111,121,120,0,0,0,0,0,0,0,0,47,111,114,114,100,56,47,53,0,0,0,0,0,0,0,0,37,115,32,119,104,105,108,101,32,111,112,101,110,105,110,103,32,37,115,10,0,0,0,0,65,108,112,104,97,0,0,0,77,99,105,114,99,108,101,0,47,111,114,114,100,56,47,52,0,0,0,0,0,0,0,0,110,111,100,101,115,101,112,0,114,111,109,97,110,0,0,0,111,117,116,112,117,116,111,114,100,101,114,0,0,0,0,0,91,32,0,0,0,0,0,0,47,111,114,114,100,56,47,51,0,0,0,0,0,0,0,0,47,111,114,114,100,56,47,50,0,0,0,0,0,0,0,0,47,111,114,114,100,56,47,49,0,0,0,0,0,0,0,0,97,103,110,101,120,116,111,117,116,101,100,103,101,0,0,0,32,37,100,32,37,100,0,0,99,109,97,112,120,58,109,97,112,0,0,0,0,0,0,0,47,111,114,114,100,55,47,55,0,0,0,0,0,0,0,0,120,100,111,116,0,0,0,0,47,111,114,114,100,55,47,54,0,0,0,0,0,0,0,0,106,112,103,58,115,118,103,0,99,105,114,99,111,0,0,0,47,111,114,114,100,55,47,53,0,0,0,0,0,0,0,0,37,108,102,44,37,108,102,0,114,97,110,107,115,101,112,0,47,98,108,117,101,115,53,47,51,0,0,0,0,0,0,0,46,46,47,46,46,47,46,46,47,108,105,98,47,110,101,97,116,111,103,101,110,47,115,109,97,114,116,95,105,110,105,95,120,46,99,0,0,0,0,0,108,101,118,101,108,115,0,0,47,111,114,114,100,55,47,52,0,0,0,0,0,0,0,0,105,100,101,97,108,95,100,105,115,116,95,115,99,104,101,109,101,32,118,97,108,117,101,32,119,114,111,110,103,0,0,0,46,46,47,46,46,47,46,46,47,112,108,117,103,105,110,47,99,111,114,101,47,103,118,114,101,110,100,101,114,95,99,111,114,101,95,118,109,108,46,99,0,0,0,0,0,0,0,0,82,73,70,70,0,0,0,0,100,105,109,103,114,97,121,0,37,37,37,37,80,97,103,101,79,114,105,101,110,116,97,116,105,111,110,58,32,37,115,10,0,0,0,0,0,0,0,0,103,114,97,121,50,53,0,0,80,97,108,97,116,105,110,111,45,66,111,108,100,0,0,0,47,111,114,114,100,55,47,51,0,0,0,0,0,0,0,0,100,101,102,97,117,108,116,100,105,115,116,0,0,0,0,0,103,105,102,58,118,109,108,0,76,111,97,100,105,110,103,32,37,115,10,0,0,0,0,0,109,100,115,0,0,0,0,0,115,97,109,112,108,101,112,111,105,110,116,115,0,0,0,0,120,121,32,112,115,101,117,100,111,45,111,114,116,104,111,103,111,110,97,108,32,99,111,110,115,116,114,97,105,110,116,115,0,0,0,0,0,0,0,0,47,111,114,114,100,55,47,50,0,0,0,0,0,0,0,0,82,79,85,78,68,40,71,68,95,98,98,40,103,41,46,76,76,46,121,41,32,61,61,32,48,0,0,0,0,0,0,0,65,103,114,97,118,101,0,0,77,115,113,117,97,114,101,0,110,111,100,101,32,37,115,44,32,112,111,115,105,116,105,111,110,32,37,115,44,32,101,120,112,101,99,116,101,100,32,116,119,111,32,100,111,117,98,108,101,115,10,0,0,0,0,0,47,111,114,114,100,55,47,49,0,0,0,0,0,0,0,0,82,76,0,0,0,0,0,0,100,111,116,116,101,100,0,0,37,108,102,44,37,108,102,44,37,108,102,44,37,108,102,44,37,108,102,0,0,0,0,0,114,101,109,111,118,101,95,111,118,101,114,108,97,112,58,32,71,114,97,112,104,118,105,122,32,110,111,116,32,98,117,105,108,116,32,119,105,116,104,32,116,114,105,97,110,103,117,108,97,116,105,111,110,32,108,105,98,114,97,114,121,10,0,0,47,111,114,114,100,54,47,54,0,0,0,0,0,0,0,0,47,111,114,114,100,54,47,53,0,0,0,0,0,0,0,0,95,76,84,88,95,108,105,98,114,97,114,121,0,0,0,0,110,111,114,109,97,108,105,122,101,0,0,0,0,0,0,0,57,58,112,114,105,115,109,0,115,112,108,105,110,101,115,32,97,110,100,32,99,108,117,115,116,101,114,32,101,100,103,101,115,32,110,111,116,32,115,117,112,112,111,114,116,101,100,32,45,32,117,115,105,110,103,32,108,105,110,101,32,115,101,103,109,101,110,116,115,10,0,0,99,95,99,110,116,32,61,61,32,48,0,0,0,0,0,0,47,111,114,114,100,54,47,52,0,0,0,0,0,0,0,0,97,103,102,105,114,115,116,111,117,116,101,100,103,101,0,0,47,111,114,114,100,54,47,51,0,0,0,0,0,0,0,0,108,101,118,101,108,32,97,115,115,105,103,110,109,101,110,116,32,99,111,110,115,116,114,97,105,110,116,115,0,0,0,0,75,80,95,82,105,103,104,116,0,0,0,0,0,0,0,0,105,110,115,116,97,108,108,95,105,110,95,114,97,110,107,44,32,108,105,110,101,32,37,100,58,32,78,68,95,111,114,100,101,114,40,37,115,41,32,91,37,100,93,32,62,32,71,68,95,114,97,110,107,40,82,111,111,116,41,91,37,100,93,46,97,110,32,91,37,100,93,10,0,0,0,0,0,0,0,0,78,68,95,110,101,120,116,40,118,41,32,61,61,32,78,85,76,76,0,0,0,0,0,0,46,46,47,46,46,47,46,46,47,108,105,98,47,100,111,116,103,101,110,47,100,111,116,105,110,105,116,46,99,0,0,0,109,97,107,101,83,112,108,105,110,101,58,32,102,97,105,108,101,100,32,116,111,32,109,97,107,101,32,115,112,108,105,110,101,32,101,100,103,101,32,40,37,115,44,37,115,41,10,0,78,68,95,105,110,40,114,105,103,104,116,41,46,115,105,122,101,32,43,32,78,68,95,111,117,116,40,114,105,103,104,116,41,46,115,105,122,101,32,61,61,32,48,0,0,0,0,0,37,115,32,45,62,32,37,115,58,32,116,97,105,108,32,105,115,32,105,110,115,105,100,101,32,104,101,97,100,32,99,108,117,115,116,101,114,32,37,115,10,0,0,0,0,0,0,0,47,111,114,114,100,54,47,50,0,0,0,0,0,0,0,0,115,45,62,115,122,32,62,32,48,0,0,0,0,0,0,0,47,111,114,114,100,54,47,49,0,0,0,0,0,0,0,0,99,97,110,110,111,116,32,102,105,110,100,32,116,114,105,97,110,103,108,101,32,112,97,116,104,0,0,0,0,0,0,0,47,98,108,117,101,115,53,47,50,0,0,0,0,0,0,0,47,111,114,114,100,53,47,53,0,0,0,0,0,0,0,0,35,37,48,50,120,37,48,50,120,37,48,50,120,0,0,0,120,109,108,0,0,0,0,0,100,101,101,112,115,107,121,98,108,117,101,0,0,0,0,0,37,37,37,37,80,97,103,101,66,111,117,110,100,105,110,103,66,111,120,58,32,37,100,32,37,100,32,37,100,32,37,100,10,0,0,0,0,0,0,0,103,114,97,121,50,48,0,0,80,97,108,97,116,105,110,111,45,82,111,109,97,110,0,0,47,111,114,114,100,53,47,52,0,0,0,0,0,0,0,0,112,110,103,58,118,109,108,0,115,104,111,114,116,112,97,116,104,0,0,0,0,0,0,0,112,111,114,116,104,111,120,121,0,0,0,0,0,0,0,0,47,111,114,114,100,53,47,51,0,0,0,0,0,0,0,0,65,99,105,114,99,0,0,0,77,100,105,97,109,111,110,100,0,0,0,0,0,0,0,0,32,114,45,120,112,32,0,0,47,111,114,114,100,53,47,50,0,0,0,0,0,0,0,0,66,84,0,0,0,0,0,0,70,65,76,83,69,0,0,0,100,97,115,104,101,100,0,0,37,108,102,44,37,108,102,44,37,108,102,44,37,91,94,44,93,37,115,0,0,0,0,0,47,111,114,114,100,53,47,49,0,0,0,0,0,0,0,0,101,112,115,102,0,0,0,0,105,110,32,114,111,117,116,101,115,112,108,105,110,101,115,44,32,105,108,108,101,103,97,108,32,118,97,108,117,101,115,32,111,102,32,112,114,101,118,32,37,100,32,97,110,100,32,110,101,120,116,32,37,100,44,32,108,105,110,101,32,37,100,10,0,0,0,0,0,0,0,0,47,111,114,114,100,52,47,52,0,0,0,0,0,0,0,0,47,111,114,114,100,52,47,51,0,0,0,0,0,0,0,0,45,45,0,0,0,0,0,0,97,103,110,101,120,116,105,110,101,100,103,101,0,0,0,0,47,111,114,114,100,52,47,50,0,0,0,0,0,0,0,0,47,111,114,114,100,52,47,49,0,0,0,0,0,0,0,0,59,10,0,0,0,0,0,0,47,111,114,114,100,51,47,51,0,0,0,0,0,0,0,0,47,98,108,117,101,115,53,47,49,0,0,0,0,0,0,0,97,114,105,97,108,0,0,0,47,111,114,114,100,51,47,50,0,0,0,0,0,0,0,0,110,111,110,101,0,0,0,0,60,63,120,109,108,0,0,0,100,101,101,112,112,105,110,107,0,0,0,0,0,0,0,0,37,37,37,37,80,97,103,101,58,32,37,100,32,37,100,10,0,0,0,0,0,0,0,0,103,114,97,121,49,53,0,0,46,46,47,46,46,47,46,46,47,108,105,98,47,99,111,109,109,111,110,47,101,109,105,116,46,99,0,0,0,0,0,0,78,101,119,67,101,110,116,117,114,121,83,99,104,108,98,107,45,66,111,108,100,73,116,97,108,105,99,0,0,0,0,0,47,111,114,114,100,51,47,49,0,0,0,0,0,0,0,0,115,117,98,115,101,116,0,0,112,111,114,116,104,111,95,121,120,0,0,0,0,0,0,0,47,111,114,97,110,103,101,115,57,47,57,0,0,0,0,0,65,97,99,117,116,101,0,0,105,110,118,104,111,117,115,101,0,0,0,0,0,0,0,0,47,111,114,97,110,103,101,115,57,47,56,0,0,0,0,0,76,82,0,0,0,0,0,0,116,114,121,105,110,103,32,116,111,32,100,101,108,101,116,101,32,97,32,110,111,110,45,108,105,110,101,10,0,0,0,0,37,108,102,44,37,108,102,44,37,108,102,44,39,37,91,94,39,93,39,0,0,0,0,0,47,111,114,97,110,103,101,115,57,47,55,0,0,0,0,0,47,111,114,97,110,103,101,115,57,47,54,0,0,0,0,0,47,111,114,97,110,103,101,115,57,47,53,0,0,0,0,0,97,103,102,105,114,115,116,105,110,101,100,103,101,0,0,0,47,111,114,97,110,103,101,115,57,47,52,0,0,0,0,0,47,111,114,97,110,103,101,115,57,47,51,0,0,0,0,0,47,111,114,97,110,103,101,115,57,47,50,0,0,0,0,0,47,98,108,117,101,115,52,47,52,0,0,0,0,0,0,0,47,111,114,97,110,103,101,115,57,47,49,0,0,0,0,0,34,32,47,62,0,0,0,0,101,112,115,0,0,0,0,0,100,97,114,107,118,105,111,108,101,116,0,0,0,0,0,0,37,37,37,37,69,110,100,80,97,103,101,58,32,37,100,10,0,0,0,0,0,0,0,0,103,114,97,121,49,48,0,0,78,101,119,67,101,110,116,117,114,121,83,99,104,108,98,107,45,82,111,109,97,110,0,0,47,111,114,97,110,103,101,115,56,47,56,0,0,0,0,0,115,118,103,58,120,100,111,116,0,0,0,0,0,0,0,0,99,105,114,99,117,105,116,0,99,111,109,112,114,101,115,115,32,37,103,32,10,0,0,0,112,115,101,117,100,111,45,111,114,116,104,111,103,111,110,97,108,32,99,111,110,115,116,114,97,105,110,116,115,0,0,0,47,111,114,97,110,103,101,115,56,47,55,0,0,0,0,0,99,111,110,106,117,103,97,116,101,95,103,114,97,100,105,101,110,116,58,32,117,110,101,120,112,101,99,116,101,100,32,108,101,110,103,116,104,32,48,32,118,101,99,116,111,114,10,0,65,69,108,105,103,0,0,0,105,110,118,116,114,97,112,101,122,105,117,109,0,0,0,0,47,111,114,97,110,103,101,115,56,47,54,0,0,0,0,0,114,97,110,107,100,105,114,0,67,101,110,116,117,114,121,32,83,99,104,111,111,108,98,111,111,107,32,76,0,0,0,0,118,105,101,119,112,111,114,116,0,0,0,0,0,0,0,0,47,111,114,97,110,103,101,115,56,47,53,0,0,0,0,0,47,111,114,97,110,103,101,115,56,47,52,0,0,0,0,0,47,111,114,97,110,103,101,115,56,47,51,0,0,0,0,0,97,103,110,101,120,116,101,100,103,101,0,0,0,0,0,0,47,111,114,97,110,103,101,115,56,47,50,0,0,0,0,0,47,111,114,97,110,103,101,115,56,47,49,0,0,0,0,0,124,101,100,103,101,108,97,98,101,108,124,0,0,0,0,0,47,111,114,97,110,103,101,115,55,47,55,0,0,0,0,0,47,98,108,117,101,115,52,47,51,0,0,0,0,0,0,0,47,111,114,97,110,103,101,115,55,47,54,0,0,0,0,0,34,32,100,97,115,104,115,116,121,108,101,61,34,100,111,116,0,0,0,0,0,0,0,0,41,10,0,0,0,0,0,0,197,208,211,198,0,0,0,0,100,97,114,107,116,117,114,113,117,111,105,115,101,0,0,0,37,37,80,97,103,101,84,114,97,105,108,101,114,10,0,0,103,114,97,121,48,53,0,0,78,101,119,67,101,110,116,117,114,121,83,99,104,108,98,107,45,73,116,97,108,105,99,0,47,111,114,97,110,103,101,115,55,47,53,0,0,0,0,0,101,112,115,58,120,100,111,116,0,0,0,0,0,0,0,0,109,111,100,101,108,0,0,0,112,111,114,116,104,111,0,0,47,111,114,97,110,103,101,115,55,47,52,0,0,0,0,0,98,122,46,115,105,122,101,0,105,110,118,116,114,105,97,110,103,108,101,0,0,0,0,0,47,111,114,97,110,103,101,115,55,47,51,0,0,0,0,0,113,117,97,110,116,117,109,0,112,97,103,101,100,105,114,61,37,115,32,105,103,110,111,114,101,100,10,0,0,0,0,0,99,111,111,114,100,115,0,0,47,111,114,97,110,103,101,115,55,47,50,0,0,0,0,0,47,111,114,97,110,103,101,115,55,47,49,0,0,0,0,0,99,99,37,115,95,37,100,0,99,111,109,112,111,117,110,100,69,100,103,101,115,58,32,99,111,117,108,100,32,110,111,116,32,99,111,110,115,116,114,117,99,116,32,111,98,115,116,97,99,108,101,115,32,45,32,102,97,108,108,105,110,103,32,98,97,99,107,32,116,111,32,115,116,114,97,105,103,104,116,32,108,105,110,101,32,101,100,103,101,115,10,0,0,0,0,0,47,111,114,97,110,103,101,115,54,47,54,0,0,0,0,0,97,103,102,105,114,115,116,101,100,103,101,0,0,0,0,0,65,103,114,97,112,104,95,116,0,0,0,0,0,0,0,0,115,97,109,101,104,101,97,100,0,0,0,0,0,0,0,0,47,111,114,97,110,103,101,115,54,47,53,0,0,0,0,0,47,111,114,97,110,103,101,115,54,47,52,0,0,0,0,0,47,111,114,97,110,103,101,115,54,47,51,0,0,0,0,0])
.concat([47,98,108,117,101,115,52,47,50,0,0,0,0,0,0,0,110,115,108,105,109,105,116,49,0,0,0,0,0,0,0,0,108,111,99,97,108,101,32,110,111,116,32,115,117,112,112,111,114,116,101,100,0,0,0,0,47,111,114,97,110,103,101,115,54,47,50,0,0,0,0,0,34,32,100,97,115,104,115,116,121,108,101,61,34,100,97,115,104,0,0,0,0,0,0,0,112,100,102,0,0,0,0,0,100,97,114,107,115,108,97,116,101,103,114,101,121,0,0,0,101,110,100,112,97,103,101,10,115,104,111,119,112,97,103,101,10,103,114,101,115,116,111,114,101,10,0,0,0,0,0,0,78,101,119,67,101,110,116,117,114,121,83,99,104,108,98,107,45,66,111,108,100,0,0,0,47,111,114,97,110,103,101,115,54,47,49,0,0,0,0,0,112,115,58,120,100,111,116,0,46,46,47,46,46,47,46,46,47,108,105,98,47,100,111,116,103,101,110,47,112,111,115,105,116,105,111,110,46,99,0,0,95,110,101,97,116,111,95,99,99,0,0,0,0,0,0,0,121,120,32,111,114,116,104,111,103,111,110,97,108,32,99,111,110,115,116,114,97,105,110,116,115,0,0,0,0,0,0,0,47,111,114,97,110,103,101,115,53,47,53,0,0,0,0,0,84,119,111,32,99,108,117,115,116,101,114,115,32,110,97,109,101,100,32,37,115,32,45,32,116,104,101,32,115,101,99,111,110,100,32,119,105,108,108,32,98,101,32,105,103,110,111,114,101,100,10,0,0,0,0,0,116,114,105,112,108,101,111,99,116,97,103,111,110,0,0,0,47,111,114,97,110,103,101,115,53,47,52,0,0,0,0,0,105,109,97,103,101,112,97,116,104,0,0,0,0,0,0,0,47,111,114,97,110,103,101,115,53,47,51,0,0,0,0,0,114,101,110,100,101,114,101,114,32,102,111,114,32,37,115,32,105,115,32,117,110,97,118,97,105,108,97,98,108,101,10,0,47,111,114,97,110,103,101,115,53,47,50,0,0,0,0,0,114,101,109,105,110,99,114,111,115,115,0,0,0,0,0,0,46,46,47,46,46,47,46,46,47,108,105,98,47,100,111,116,103,101,110,47,102,108,97,116,46,99,0,0,0,0,0,0,47,111,114,97,110,103,101,115,53,47,49,0,0,0,0,0,97,103,110,120,116,111,117,116,0,0,0,0,0,0,0,0,47,111,114,97,110,103,101,115,52,47,52,0,0,0,0,0,47,111,114,97,110,103,101,115,52,47,51,0,0,0,0,0,46,46,47,46,46,47,46,46,47,108,105,98,47,100,111,116,103,101,110,47,102,97,115,116,103,114,46,99,0,0,0,0,47,111,114,97,110,103,101,115,52,47,50,0,0,0,0,0,47,98,108,117,101,115,52,47,49,0,0,0,0,0,0,0,47,111,114,97,110,103,101,115,52,47,49,0,0,0,0,0,34,32,119,101,105,103,104,116,61,34,37,46,48,102,112,116,0,0,0,0,0,0,0,0,37,80,68,70,45,0,0,0,100,97,114,107,115,108,97,116,101,103,114,97,121,0,0,0,48,32,48,32,48,32,101,100,103,101,99,111,108,111,114,10,0,0,0,0,0,0,0,0,66,111,111,107,109,97,110,45,68,101,109,105,73,116,97,108,105,99,0,0,0,0,0,0,47,111,114,97,110,103,101,115,51,47,51,0,0,0,0,0,106,112,103,58,120,100,111,116,0,0,0,0,0,0,0,0,37,115,32,97,116,116,114,105,98,117,116,101,32,118,97,108,117,101,32,109,117,115,116,32,98,101,32,49,32,111,114,32,50,32,45,32,105,103,110,111,114,105,110,103,10,0,0,0,101,100,103,101,32,108,97,98,101,108,115,32,119,105,116,104,32,115,112,108,105,110,101,115,61,99,117,114,118,101,100,32,110,111,116,32,115,117,112,112,111,114,116,101,100,32,105,110,32,100,111,116,32,45,32,117,115,101,32,120,108,97,98,101,108,115,10,0,0,0,0,0,111,114,116,104,111,121,120,0,47,111,114,97,110,103,101,115,51,47,50,0,0,0,0,0,100,111,117,98,108,101,111,99,116,97,103,111,110,0,0,0,47,111,114,97,110,103,101,115,51,47,49,0,0,0,0,0,71,68,70,79,78,84,80,65,84,72,61,0,0,0,0,0,112,104,97,115,101,0,0,0,108,97,121,111,117,116,32,119,97,115,32,110,111,116,32,100,111,110,101,10,0,0,0,0,99,111,110,99,101,110,116,114,97,116,101,61,116,114,117,101,32,109,97,121,32,110,111,116,32,119,111,114,107,32,99,111,114,114,101,99,116,108,121,46,10,0,0,0,0,0,0,0,47,103,114,101,121,115,57,47,57,0,0,0,0,0,0,0,47,103,114,101,121,115,57,47,56,0,0,0,0,0,0,0,108,104,101,97,100,0,0,0,47,103,114,101,121,115,57,47,55,0,0,0,0,0,0,0,97,103,102,115,116,111,117,116,0,0,0,0,0,0,0,0,37,115,32,119,97,115,32,97,108,114,101,97,100,121,32,105,110,32,97,32,114,97,110,107,115,101,116,44,32,100,101,108,101,116,101,100,32,102,114,111,109,32,99,108,117,115,116,101,114,32,37,115,10,0,0,0,46,46,47,46,46,47,46,46,47,108,105,98,47,100,111,116,103,101,110,47,99,108,97,115,115,50,46,99,0,0,0,0,47,103,114,101,121,115,57,47,54,0,0,0,0,0,0,0,47,103,114,101,121,115,57,47,53,0,0,0,0,0,0,0,47,103,114,101,121,115,57,47,52,0,0,0,0,0,0,0,46,46,47,46,46,47,46,46,47,108,105,98,47,99,105,114,99,111,103,101,110,47,110,111,100,101,108,105,115,116,46,99,0,0,0,0,0,0,0,0,47,103,114,101,121,115,57,47,51,0,0,0,0,0,0,0,47,98,108,117,101,115,51,47,51,0,0,0,0,0,0,0,60,118,58,115,116,114,111,107,101,32,99,111,108,111,114,61,34,0,0,0,0,0,0,0,35,32,71,101,110,101,114,97,116,101,100,32,98,121,32,0,106,112,101,103,0,0,0,0,100,97,114,107,115,108,97,116,101,98,108,117,101,0,0,0,37,37,32,37,115,10,0,0,66,111,111,107,109,97,110,45,76,105,103,104,116,0,0,0,47,103,114,101,121,115,57,47,50,0,0,0,0,0,0,0,32,50,10,0,0,0,0,0,106,112,101,58,120,100,111,116,0,0,0,0,0,0,0,0,46,46,47,46,46,47,46,46,47,108,105,98,47,99,105,114,99,111,103,101,110,47,100,101,103,108,105,115,116,46,99,0,115,116,114,101,115,115,119,116,0,0,0,0,0,0,0,0,120,121,32,111,114,116,104,111,103,111,110,97,108,32,99,111,110,115,116,114,97,105,110,116,115,0,0,0,0,0,0,0,47,103,114,101,121,115,57,47,49,0,0,0,0,0,0,0,111,110,101,98,108,111,99,107,0,0,0,0,0,0,0,0,100,111,117,98,108,101,99,105,114,99,108,101,0,0,0,0,47,103,114,101,121,115,56,47,56,0,0,0,0,0,0,0,68,79,84,70,79,78,84,80,65,84,72,0,0,0,0,0,37,0,0,0,73,0,0,0,58,0,0,0,37,0,0,0,77,0,0,0,58,0,0,0,37,0,0,0,83,0,0,0,32,0,0,0,37,0,0,0,112,0,0,0,0,0,0,0,95,115,112,97,110,95,37,100,0,0,0,0,0,0,0,0,47,103,114,101,121,115,56,47,55,0,0,0,0,0,0,0,47,103,114,101,121,115,56,47,54,0,0,0,0,0,0,0,47,103,114,101,121,115,56,47,53,0,0,0,0,0,0,0,97,103,110,120,116,105,110,0,47,103,114,101,121,115,56,47,52,0,0,0,0,0,0,0,103,114,97,112,104,118,105,122,0,0,0,0,0,0,0,0,108,105,98,112,97,116,104,47,37,115,58,37,100,58,32,37,115,10,0,0,0,0,0,0,47,103,114,101,121,115,56,47,51,0,0,0,0,0,0,0,47,103,114,101,121,115,56,47,50,0,0,0,0,0,0,0,46,46,47,46,46,47,46,46,47,108,105,98,47,112,97,116,104,112,108,97,110,47,99,118,116,46,99,0,0,0,0,0,47,98,108,117,101,115,51,47,50,0,0,0,0,0,0,0,37,73,58,37,77,58,37,83,32,37,112,0,0,0,0,0,47,103,114,101,121,115,56,47,49,0,0,0,0,0,0,0,60,47,118,58,115,104,97,112,101,62,10,0,0,0,0,0,32,80,97,103,101,115,58,32,37,100,10,0,0,0,0,0,255,216,255,224,0,0,0,0,100,97,114,107,115,101,97,103,114,101,101,110,0,0,0,0,103,115,97,118,101,10,0,0,102,108,101,115,104,0,0,0,66,111,111,107,109,97,110,45,76,105,103,104,116,73,116,97,108,105,99,0,0,0,0,0,47,103,114,101,121,115,55,47,55,0,0,0,0,0,0,0,34,32,110,97,109,101,61,34,0,0,0,0,0,0,0,0,49,50,48,48,0,0,0,0,106,112,101,103,58,120,100,111,116,0,0,0,0,0,0,0,108,97,121,111,117,116,32,97,98,111,114,116,101,100,10,0,111,114,116,104,111,120,121,0,47,103,114,101,121,115,55,47,54,0,0,0,0,0,0,0,47,97,99,99,101,110,116,51,47,50,0,0,0,0,0,0,115,113,117,97,114,101,0,0,47,103,114,101,121,115,55,47,53,0,0,0,0,0,0,0,102,111,110,116,112,97,116,104,0,0,0,0,0,0,0,0,76,97,121,111,117,116,32,119,97,115,32,110,111,116,32,100,111,110,101,46,32,32,77,105,115,115,105,110,103,32,108,97,121,111,117,116,32,112,108,117,103,105,110,115,63,32,10,0,115,116,100,58,58,98,97,100,95,97,108,108,111,99,0,0,47,103,114,101,121,115,55,47,52,0,0,0,0,0,0,0,47,103,114,101,121,115,55,47,51,0,0,0,0,0,0,0,47,103,114,101,121,115,55,47,50,0,0,0,0,0,0,0,97,103,102,115,116,105,110,0,115,105,103,110,101,100,32,99,104,97,114,0,0,0,0,0,47,103,114,101,121,115,55,47,49,0,0,0,0,0,0,0,47,103,114,101,121,115,54,47,54,0,0,0,0,0,0,0,38,35,52,53,59,0,0,0,47,103,114,101,121,115,54,47,53,0,0,0,0,0,0,0,47,98,108,117,101,115,51,47,49,0,0,0,0,0,0,0,37,0,0,0,97,0,0,0,32,0,0,0,37,0,0,0,98,0,0,0,32,0,0,0,37,0,0,0,100,0,0,0,32,0,0,0,37,0,0,0,72,0,0,0,58,0,0,0,37,0,0,0,77,0,0,0,58,0,0,0,37,0,0,0,83,0,0,0,32,0,0,0,37,0,0,0,89,0,0,0,0,0,0,0,0,0,0,0,47,103,114,101,121,115,54,47,52,0,0,0,0,0,0,0,34,47,62,0,0,0,0,0,32,84,105,116,108,101,58,32,0,0,0,0,0,0,0,0,103,105,102,0,0,0,0,0,100,97,114,107,115,97,108,109,111,110,0,0,0,0,0,0,103,114,101,115,116,111,114,101,10,0,0,0,0,0,0,0,66,111,111,107,109,97,110,45,68,101,109,105,0,0,0,0,47,103,114,101,121,115,54,47,51,0,0,0,0,0,0,0,60,109,97,112,32,105,100,61,34,0,0,0,0,0,0,0,45,50,10,0,0,0,0,0,103,105,102,58,120,100,111,116,0,0,0,0,0,0,0,0,97,110,116,105,113,117,101,119,104,105,116,101,0,0,0,0,108,101,118,101,108,115,103,97,112,0,0,0,0,0,0,0,111,114,116,104,111,95,121,120,0,0,0,0,0,0,0,0,47,103,114,101,121,115,54,47,50,0,0,0,0,0,0,0,85,110,107,110,111,119,110,32,34,115,112,108,105,110,101,115,34,32,118,97,108,117,101,58,32,34,37,115,34,32,45,32,105,103,110,111,114,101,100,10,0,0,0,0,0,0,0,0,114,101,99,116,97,110,103,108,101,0,0,0,0,0,0,0,47,103,114,101,121,115,54,47,49,0,0,0,0,0,0,0,99,111,110,100,101,110,115,101,100,0,0,0,0,0,0,0,101,112,115,58,112,115,0,0,47,103,114,101,121,115,53,47,53,0,0,0,0,0,0,0,47,103,114,101,121,115,53,47,52,0,0,0,0,0,0,0,112,111,118,58,112,111,118,0,47,103,114,101,121,115,53,47,51,0,0,0,0,0,0,0,97,103,110,120,116,101,100,103,101,0,0,0,0,0,0,0,46,92,34,32,0,0,0,0,105,109,97,112,58,109,97,112,0,0,0,0,0,0,0,0,47,103,114,101,121,115,53,47,50,0,0,0,0,0,0,0,100,111,116,0,0,0,0,0,47,103,114,101,121,115,53,47,49,0,0,0,0,0,0,0,106,112,101,58,115,118,103,0,116,119,111,112,105,0,0,0,65,103,110,111,100,101,105,110,102,111,95,116,0,0,0,0,37,46,48,51,102,0,0,0,47,103,114,101,121,115,52,47,52,0,0,0,0,0,0,0,115,99,97,108,101,0,0,0,47,97,99,99,101,110,116,56,47,56,0,0,0,0,0,0,114,101,112,117,108,115,105,118,101,102,111,114,99,101,0,0,37,97,32,37,98,32,37,100,32,37,72,58,37,77,58,37,83,32,37,89,0,0,0,0,47,103,114,101,121,115,52,47,51,0,0,0,0,0,0,0,110,122,32,62,32,48,0,0,32,101,32,0,0,0,0,0,35,0,0,0,0,0,0,0,71,73,70,56,0,0,0,0,100,97,114,107,114,101,100,0,32,32,47,66,111,114,100,101,114,32,91,32,48,32,48,32,48,32,93,10,32,32,47,65,99,116,105,111,110,32,60,60,32,47,83,117,98,116,121,112,101,32,47,85,82,73,32,47,85,82,73,32,37,115,32,62,62,10,32,32,47,83,117,98,116,121,112,101,32,47,76,105,110,107,10,47,65,78,78,32,112,100,102,109,97,114,107,10,0,0,0,0,0,0,0,0,102,101,108,100,115,112,97,114,0,0,0,0,0,0,0,0,47,103,114,101,121,115,52,47,50,0,0,0,0,0,0,0,84,105,109,101,115,45,73,116,97,108,105,99,0,0,0,0,68,97,109,112,105,110,103,0,83,105,110,103,108,101,10,0,112,110,103,58,120,100,111,116,0,0,0,0,0,0,0,0,105,115,32,117,110,100,101,102,105,110,101,100,46,32,82,101,118,101,114,116,105,110,103,32,116,111,32,116,104,101,32,115,104,111,114,116,101,115,116,32,112,97,116,104,32,109,111,100,101,108,46,10,0,0,0,0,67,111,117,108,100,32,110,111,116,32,108,111,97,100,32,34,37,115,34,32,45,32,37,115,10,0,0,0,0,0,0,0,109,97,107,101,80,111,108,121,58,32,117,110,107,110,111,119,110,32,115,104,97,112,101,32,116,121,112,101,32,37,115,10,0,0,0,0,0,0,0,0,111,114,116,104,111,103,111,110,97,108,32,99,111,110,115,116,114,97,105,110,116,115,0,0,47,103,114,101,121,115,52,47,49,0,0,0,0,0,0,0,82,79,85,78,68,40,71,68,95,98,98,40,103,41,46,76,76,46,120,41,32,61,61,32,48,0,0,0,0,0,0,0,101,115,0,0,0,0,0,0,114,101,99,116,0,0,0,0,37,46,53,103,44,37,46,53,103,44,37,46,53,103,44,37,46,53,103,32,0,0,0,0,47,103,114,101,121,115,51,47,51,0,0,0,0,0,0,0,117,110,109,97,116,99,104,101,100,32,39,40,39,32,105,110,32,115,116,121,108,101,58,32,37,115,10,0,0,0,0,0,47,103,114,101,121,115,51,47,50,0,0,0,0,0,0,0,99,103,0,0,0,0,0,0,47,103,114,101,121,115,51,47,49,0,0,0,0,0,0,0,102,100,112,32,100,111,101,115,32,110,111,116,32,115,117,112,112,111,114,116,32,115,116,97,114,116,61,115,101,108,102,32,45,32,105,103,110,111,114,105,110,103,10,0,0,0,0,0,46,46,47,46,46,47,46,46,47,108,105,98,47,102,100,112,103,101,110,47,99,111,109,112,46,99,0,0,0,0,0,0,47,103,114,101,101,110,115,57,47,57,0,0,0,0,0,0,97,103,102,115,116,101,100,103,101,0,0,0,0,0,0,0,82,105,103,104,116,0,0,0,71,68,95,114,97,110,107,40,103,41,91,114,93,46,110,32,60,61,32,71,68,95,114,97,110,107,40,103,41,91,114,93,46,97,110,0,0,0,0,0,47,103,114,101,101,110,115,57,47,56,0,0,0,0,0,0,120,120,120,0,0,0,0,0,46,46,47,46,46,47,46,46,47,108,105,98,47,100,111,116,103,101,110,47,99,111,110,99,46,99,0,0,0,0,0,0,37,115,32,45,62,32,37,115,58,32,104,101,97,100,32,110,111,116,32,105,110,115,105,100,101,32,104,101,97,100,32,99,108,117,115,116,101,114,32,37,115,10,0,0,0,0,0,0,47,103,114,101,101,110,115,57,47,55,0,0,0,0,0,0,69,68,95,116,111,95,118,105,114,116,40,111,114,105,103,41,32,33,61,32,78,85,76,76,0,0,0,0,0,0,0,0,103,118,119,114,105,116,101,95,110,111,95,122,32,112,114,111,98,108,101,109,32,37,100,10,0,0,0,0,0,0,0,0,97,114,116,105,99,117,108,97,116,105,111,110,95,112,111,115,0,0,0,0,0,0,0,0,46,46,47,46,46,47,46,46,47,108,105,98,47,99,105,114,99,111,103,101,110,47,98,108,111,99,107,116,114,101,101,46,99,0,0,0,0,0,0,0,95,99,108,111,110,101,95,37,100,0,0,0,0,0,0,0,47,103,114,101,101,110,115,57,47,54,0,0,0,0,0,0,100,101,115,116,105,110,97,116,105,111,110,32,112,111,105,110,116,32,110,111,116,32,105,110,32,97,110,121,32,116,114,105,97,110,103,108,101,0,0,0,47,97,99,99,101,110,116,56,47,55,0,0,0,0,0,0,47,103,114,101,101,110,115,57,47,53,0,0,0,0,0,0,32,108,32,0,0,0,0,0,32,45,97,110,99,104,111,114,32,101,0,0,0,0,0,0,98,109,112,0,0,0,0,0,100,97,114,107,111,114,99,104,105,100,0,0,0,0,0,0,32,93,10,0,0,0,0,0,100,117,115,116,121,114,111,115,101,0,0,0,0,0,0,0,47,103,114,101,101,110,115,57,47,52,0,0,0,0,0,0,72,101,108,118,101,116,105,99,97,45,78,97,114,114,111,119,45,66,111,108,100,79,98,108,105,113,117,101,0,0,0,0,100,101,102,97,117,108,116,32,0,0,0,0,0,0,0,0,49,48,48,46,48,48,10,0,115,118,103,58,100,111,116,0,111,114,116,104,111,0,0,0,47,103,114,101,101,110,115,57,47,51,0,0,0,0,0,0,114,117,101,0,0,0,0,0,82,0,0,0,0,0,0,0,99,111,109,112,111,110,101,110,116,0,0,0,0,0,0,0,37,46,53,103,44,37,46,53,103,44,37,46,53,103,44,37,46,53,103,0,0,0,0,0,47,103,114,101,101,110,115,57,47,50,0,0,0,0,0,0,103,100,0,0,0,0,0,0,103,101,116,115,112,108,105,110,101,112,111,105,110,116,115,58,32,110,111,32,115,112,108,105,110,101,32,112,111,105,110,116,115,32,97,118,97,105,108,97,98,108,101,32,102,111,114,32,101,100,103,101,32,40,37,115,44,37,115,41,10,0,0,0,37,0,0,0,72,0,0,0,58,0,0,0,37,0,0,0,77,0,0,0,58,0,0,0,37,0,0,0,83,0,0,0,0,0,0,0,0,0,0,0,116,114,117,110,99,97,116,105,110,103,32,115,116,121,108,101,32,39,37,115,39,10,0,0,47,103,114,101,101,110,115,57,47,49,0,0,0,0,0,0,115,104,97,112,101,102,105,108,101,0,0,0,0,0,0,0,105,110,32,114,111,117,116,101,115,112,108,105,110,101,115,44,32,99,97,110,110,111,116,32,102,105,110,100,32,78,79,82,77,65,76,32,101,100,103,101,10,0,0,0,0,0,0,0,99,97,110,39,116,32,102,105,110,100,32,108,105,98,114,97,114,121,32,102,105,108,101,32,37,115,10,0,0,0,0,0,102,111,114,99,101,108,97,98,101,108,115,0,0,0,0,0,37,100,32,37,100,32,37,100,32,37,100,0,0,0,0,0,47,103,114,101,101,110,115,56,47,56,0,0,0,0,0,0,115,111,108,105,100,0,0,0,46,46,47,46,46,47,46,46,47,108,105,98,47,99,111,109,109,111,110,47,117,116,105,108,115,46,99,0,0,0,0,0,45,62,0,0,0,0,0,0,47,103,114,101,101,110,115,56,47,55,0,0,0,0,0,0,97,103,101,100,103,101,0,0,47,103,114,101,101,110,115,56,47,54,0,0,0,0,0,0,47,103,114,101,101,110,115,56,47,53,0,0,0,0,0,0,46,46,47,46,46,47,46,46,47,108,105,98,47,99,111,109,109,111,110,47,115,112,108,105,110,101,115,46,99,0,0,0,71,86,66,73,78,68,73,82,0,0,0,0,0,0,0,0,47,103,114,101,101,110,115,56,47,52,0,0,0,0,0,0,47,97,99,99,101,110,116,56,47,54,0,0,0,0,0,0,47,103,114,101,101,110,115,56,47,51,0,0,0,0,0,0,37,46,48,102,44,37,46,48,102,32,0,0,0,0,0,0,32,45,97,110,99,104,111,114,32,119,0,0,0,0,0,0,66,77,0,0,0,0,0,0,100,97,114,107,111,114,97,110,103,101,0,0,0,0,0,0,91,32,47,82,101,99,116,32,91,32,0,0,0,0,0,0,100,107,103,114,101,101,110,99,111,112,112,101,114,0,0,0,110,111,32,109,101,109,111,114,121,32,102,114,111,109,32,122,109,97,108,108,111,99,40,41,10,0,0,0,0,0,0,0,47,103,114,101,101,110,115,56,47,50,0,0,0,0,0,0,72,101,108,118,101,116,105,99,97,45,78,97,114,114,111,119,0,0,0,0,0,0,0,0,98,97,115,101,32,114,101,102,101,114,101,114,10,0,0,0,76,101,116,116,101,114,10,0,37,115,37,100,32,45,0,0,101,112,115,58,100,111,116,0,120,32,97,110,100,32,121,32,115,99,97,108,105,110,103,0,47,103,114,101,101,110,115,56,47,49,0,0,0,0,0,0,112,108,105,110,101,0,0,0,81,0,0,0,0,0,0,0,98,111,120,51,100,0,0,0,47,103,114,101,101,110,115,55,47,55,0,0,0,0,0,0,111,117,116,32,111,102,32,100,121,110,97,109,105,99,32,109,101,109,111,114,121,32,105,110,32,97,97,103,95,99,114,101,97,116,101,95,98,117,102,102,101,114,40,41,0,0,0,0,117,110,109,97,116,99,104,101,100,32,39,41,39,32,105,110,32,115,116,121,108,101,58,32,37,115,10,0,0,0,0,0,47,103,114,101,101,110,115,55,47,54,0,0,0,0,0,0,47,103,114,101,101,110,115,55,47,53,0,0,0,0,0,0,47,103,114,101,101,110,115,55,47,52,0,0,0,0,0,0,97,103,100,101,103,114,101,101,0,0,0,0,0,0,0,0,47,103,114,101,101,110,115,55,47,51,0,0,0,0,0,0,97,103,114,97,112,104,111,102,32,97,32,98,97,100,32,111,98,106,101,99,116,0,0,0,47,103,114,101,101,110,115,55,47,50,0,0,0,0,0,0,47,103,114,101,101,110,115,55,47,49,0,0,0,0,0,0,47,97,99,99,101,110,116,56,47,53,0,0,0,0,0,0,47,103,114,101,101,110,115,54,47,54,0,0,0,0,0,0,32,109,32,0,0,0,0,0,32,37,100,125,0,0,0,0,100,97,114,107,111,108,105,118,101,103,114,101,101,110,0,0,32,37,115,32,97,108,105,103,110,101,100,116,101,120,116,10,0,0,0,0,0,0,0,0,100,97,114,107,119,111,111,100,0,0,0,0,0,0,0,0,47,103,114,101,101,110,115,54,47,53,0,0,0,0,0,0,72,101,108,118,101,116,105,99,97,45,78,97,114,114,111,119,45,79,98,108,105,113,117,101,0,0,0,0,0,0,0,0,60,47,109,97,112,62,10,0,73,110,99,104,101,115,10,0,95,116,108,100,114,97,119,95,0,0,0,0,0,0,0,0,112,115,58,100,111,116,0,0,114,111,117,116,101,115,112,108,105,110,101,115,105,110,105,116,58,32,99,97,110,110,111,116,32,97,108,108,111,99,97,116,101,32,112,115,10,0,0,0,115,99,97,108,105,110,103,32,102,97,99,116,111,114,32,61,32,37,102,10,0,0,0,0,100,105,115,116,32,62,32,48,0,0,0,0,0,0,0,0,115,99,97,108,101,120,121,0,47,103,114,101,101,110,115,54,47,52,0,0,0,0,0,0,111,108,121,108,105,110,101,0,102,111,108,100,101,114,0,0,101,44,37,46,53,103,44,37,46,53,103,32,0,0,0,0,47,103,114,101,101,110,115,54,47,51,0,0,0,0,0,0,103,108,111,98,97,108,0,0,37,0,0,0,109,0,0,0,47,0,0,0,37,0,0,0,100,0,0,0,47,0,0,0,37,0,0,0,121,0,0,0,0,0,0,0,0,0,0,0,110,101,115,116,105,110,103,32,110,111,116,32,97,108,108,111,119,101,100,32,105,110,32,115,116,121,108,101,58,32,37,115,10,0,0,0,0,0,0,0,47,103,114,101,101,110,115,54,47,50,0,0,0,0,0,0,47,103,114,101,101,110,115,54,47,49,0,0,0,0,0,0,47,103,114,101,101,110,115,53,47,53,0,0,0,0,0,0,97,103,100,101,108,110,111,100,101,0,0,0,0,0,0,0,47,103,114,101,101,110,115,53,47,52,0,0,0,0,0,0,47,103,114,101,101,110,115,53,47,51,0,0,0,0,0,0,47,103,114,101,101,110,115,53,47,50,0,0,0,0,0,0,103,114,97,112,104,32,0,0,47,97,99,99,101,110,116,56,47,52,0,0,0,0,0,0,58,0,0,0,0,0,0,0,47,103,114,101,101,110,115,53,47,49,0,0,0,0,0,0,60,118,58,112,97,116,104,32,118,61,34,0,0,0,0,0,37,33,80,83,45,65,100,111,98,101,45,0,0,0,0,0,100,97,114,107,109,97,103,101,110,116,97,0,0,0,0,0,37,115,37,115,37,115,0,0,32,109,111,118,101,116,111,32,0,0,0,0,0,0,0,0,47,103,114,101,101,110,115,52,47,52,0,0,0,0,0,0,72,101,108,118,101,116,105,99,97,45,78,97,114,114,111,119,45,66,111,108,100,0,0,0,34,62,10,0,0,0,0,0,67,101,110,116,101,114,10,0,95,104,108,100,114,97,119,95,0,0,0,0,0,0,0,0,106,112,103,58,100,111,116,0,37,100,32,37,100,10,0,0,113,116,49,45,62,110,32,62,32,48,32,38,38,32,113,116,50,45,62,110,32,62,32,48,0,0,0,0,0,0,0,0,111,117,116,32,111,102,32,109,101,109,111,114,121,10,0,0,102,32,60,32,103,114,97,112,104,91,106,93,46,110,101,100,103,101,115,0,0,0,0,0,111,108,100,32,115,99,97,108,105,110,103,0,0,0,0,0,47,103,114,101,101,110,115,52,47,51,0,0,0,0,0,0,114,116,104,111,0,0,0,0,115,44,37,46,53,103,44,37,46,53,103,32,0,0,0,0,47,103,114,101,101,110,115,52,47,50,0,0,0,0,0,0,108,111,99,97,108,0,0,0,67,111,117,114,105,101,114,45,79,98,108,105,113,117,101,0,105,110,32,99,108,117,115,116,101,114,32,37,115,10,0,0,47,103,114,101,101,110,115,52,47,49,0,0,0,0,0,0,115,116,121,108,101,0,0,0,105,110,32,108,97,98,101,108,32,111,102,32,103,114,97,112,104,32,37,115,10,0,0,0,47,103,114,101,101,110,115,51,47,51,0,0,0,0,0,0,47,103,114,101,101,110,115,51,47,50,0,0,0,0,0,0,97,103,108,97,115,116,110,111,100,101,0,0,0,0,0,0,71,86,67,95,116,0,0,0,47,103,114,101,101,110,115,51,47,49,0,0,0,0,0,0,104,101,97,100,32,110,111,100,101,32,37,115,32,105,110,115,105,100,101,32,116,97,105,108,32,99,108,117,115,116,101,114,32,37,115,10,0,0,0,0,47,103,110,98,117,57,47,57,0,0,0,0,0,0,0,0,116,97,105,108,32,110,111,100,101,32,37,115,32,105,110,115,105,100,101,32,104,101,97,100,32,99,108,117,115,116,101,114,32,37,115,10,0,0,0,0,47,103,110,98,117,57,47,56,0,0,0,0,0,0,0,0,47,97,99,99,101,110,116,56,47,51,0,0,0,0,0,0,104,101,97,100,32,99,108,117,115,116,101,114,32,37,115,32,105,110,115,105,100,101,32,116,97,105,108,32,99,108,117,115,116,101,114,32,37,115,10,0,47,103,110,98,117,57,47,55,0,0,0,0,0,0,0,0,32,119,105,100,116,104,58,32,37,100,59,32,104,101,105,103,104,116,58,32,37,100,34,32,102,105,108,108,101,100,61,34,102,97,108,115,101,34,62,0,32,45,102,111,110,116,32,123,0,0,0,0,0,0,0,0,112,110,103,0,0,0,0,0,100,97,114,107,107,104,97,107,105,0,0,0,0,0,0,0,116,97,105,108,32,99,108,117,115,116,101,114,32,37,115,32,105,110,115,105,100,101,32,104,101,97,100,32,99,108,117,115,116,101,114,32,37,115,10,0,32,47,37,115,32,115,101,116,95,102,111,110,116,10,0,0,100,97,114,107,116,97,110,0,47,103,110,98,117,57,47,54,0,0,0,0,0,0,0,0,72,101,108,118,101,116,105,99,97,45,66,111,108,100,79,98,108,105,113,117,101,0,0,0,80,111,114,116,114,97,105,116,10,0,0,0,0,0,0,0,95,116,100,114,97,119,95,0,106,112,101,58,100,111,116,0,120,120,120,32,37,100,32,37,100,10,0,0,0,0,0,0,113,116,50,45,62,110,32,62,32,48,0,0,0,0,0,0,99,108,117,115,116,101,114,32,99,121,99,108,101,32,37,115,32,45,45,32,37,115,32,110,111,116,32,115,117,112,112,111,114,116,101,100,10,0,0,0,100,105,103,114,97,112,104,0,111,115,99,97,108,101,0,0,47,103,110,98,117,57,47,53,0,0,0,0,0,0,0,0,110,97,109,101,0,0,0,0,110,111,116,101,0,0,0,0,37,46,53,103,32,37,46,53,103,0,0,0,0,0,0,0,47,103,110,98,117,57,47,52,0,0,0,0,0,0,0,0,122,119,110,106,0,0,0,0,47,103,110,98,117,57,47,51,0,0,0,0,0,0,0,0,122,119,106,0,0,0,0,0,47,103,110,98,117,57,47,50,0,0,0,0,0,0,0,0,122,101,116,97,0,0,0,0,47,103,110,98,117,57,47,49,0,0,0,0,0,0,0,0,97,103,112,114,101,118,110,111,100,101,0,0,0,0,0,0,121,117,109,108,0,0,0,0,47,103,110,98,117,56,47,56,0,0,0,0,0,0,0,0,47,103,110,98,117,56,47,55,0,0,0,0,0,0,0,0,121,101,110,0,0,0,0,0,121,97,99,117,116,101,0,0,47,103,110,98,117,56,47,54,0,0,0,0,0,0,0,0,47,97,99,99,101,110,116,56,47,50,0,0,0,0,0,0,120,105,0,0,0,0,0,0,47,103,110,98,117,56,47,53,0,0,0,0,0,0,0,0,32,60,118,58,115,104,97,112,101,32,115,116,121,108,101,61,34,112,111,115,105,116,105,111,110,58,97,98,115,111,108,117,116,101,59,32,0,0,0,0,137,80,78,71,13,10,26,10,0,0,0,0,0,0,0,0,100,97,114,107,103,114,101,121,0,0,0,0,0,0,0,0,119,101,105,101,114,112,0,0,32,101,108,108,105,112,115,101,95,112,97,116,104,32,115,116,114,111,107,101,10,0,0,0,47,103,110,98,117,56,47,52,0,0,0,0,0,0,0,0,72,101,108,118,101,116,105,99,97,45,79,98,108,105,113,117,101,0,0,0,0,0,0,0,44,37,100,44,37,100,0,0,35,32,80,97,103,101,115,58,32,37,100,10,0,0,0,0,95,104,100,114,97,119,95,0,106,112,101,103,58,100,111,116,0,0,0,0,0,0,0,0,119,103,116,32,62,32,48,0,117,117,109,108,0,0,0,0,100,105,114,0,0,0,0,0,105,112,115,101,112,0,0,0,47,103,110,98,117,56,47,51,0,0,0,0,0,0,0,0,117,112,115,105,108,111,110,0,111,99,116,97,103,111,110,0,47,103,110,98,117,56,47,50,0,0,0,0,0,0,0,0,102,0,0,0,97,0,0,0,108,0,0,0,115,0,0,0,101,0,0,0,0,0,0,0,117,112,115,105,104,0,0,0,47,103,110,98,117,56,47,49,0,0,0,0,0,0,0,0,117,109,108,0,0,0,0,0,47,103,110,98,117,55,47,55,0,0,0,0,0,0,0,0,117,103,114,97,118,101,0,0,47,103,110,98,117,55,47,54,0,0,0,0,0,0,0,0,97,103,110,101,120,116,110,111,100,101,0,0,0,0,0,0,117,99,105,114,99,0,0,0,47,103,110,98,117,55,47,53,0,0,0,0,0,0,0,0,117,97,114,114,0,0,0,0,47,103,110,98,117,55,47,52,0,0,0,0,0,0,0,0,117,97,99,117,116,101,0,0,47,103,110,98,117,55,47,51,0,0,0,0,0,0,0,0,47,97,99,99,101,110,116,56,47,49,0,0,0,0,0,0,117,65,114,114,0,0,0,0,47,103,110,98,117,55,47,50,0,0,0,0,0,0,0,0,32,45,45,62,10,0,0,0,32,45,116,101,120,116,32,123,0,0,0,0,0,0,0,0,47,103,110,98,117,55,47,49,0,0,0,0,0,0,0,0,40,108,105,98,41,0,0,0,100,97,114,107,103,114,101,101,110,0,0,0,0,0,0,0,116,114,97,100,101,0,0,0,32,101,108,108,105,112,115,101,95,112,97,116,104,32,102,105,108,108,10,0,0,0,0,0,72,101,108,118,101,116,105,99,97,45,66,111,108,100,0,0,37,100,44,37,100,0,0,0,35,32,84,105,116,108,101,58,32,37,115,10,0,0,0,0,95,108,100,114,97,119,95,0,103,105,102,58,100,111,116,0,113,45,62,108,0,0,0,0,116,105,109,101,115,0,0,0,78,68,95,105,100,40,110,112,41,32,61,61,32,105,0,0,115,117,98,103,114,97,112,104,0,0,0,0,0,0,0,0,118,112,115,99,0,0,0,0,47,103,110,98,117,54,47,54,0,0,0,0,0,0,0,0,105,110,101,0,0,0,0,0,115,101,112,116,97,103,111,110,0,0,0,0,0,0,0,0,47,103,110,98,117,54,47,53,0,0,0,0,0,0,0,0,109,111,110,111,115,112,97,99,101,0,0,0,0,0,0,0,102,97,108,115,101,0,0,0,98,103,99,111,108,111,114,0,116,104,111,114,110,0,0,0,47,103,110,98,117,54,47,52,0,0,0,0,0,0,0,0,116,104,105,110,115,112,0,0,47,103,110,98,117,54,47,51,0,0,0,0,0,0,0,0,116,104,101,116,97,115,121,109,0,0,0,0,0,0,0,0,47,103,110,98,117,54,47,50,0,0,0,0,0,0,0,0,97,103,102,105,114,115,116,110,111,100,101,0,0,0,0,0,116,104,101,116,97,0,0,0,47,103,110,98,117,54,47,49,0,0,0,0,0,0,0,0,116,104,101,114,101,52,0,0,47,103,110,98,117,53,47,53,0,0,0,0,0,0,0,0,116,97,117,0,0,0,0,0,47,103,110,98,117,53,47,52,0,0,0,0,0,0,0,0,47,97,99,99,101,110,116,55,47,55,0,0,0,0,0,0,115,122,108,105,103,0,0,0,47,103,110,98,117,53,47,51,0,0,0,0,0,0,0,0,32,32,32,32,32,32,60,33,45,45,32,0,0,0,0,0,32,99,114,101,97,116,101,32,116,101,120,116,32,0,0,0,119,101,98,112,0,0,0,0,100,97,114,107,103,114,97,121,0,0,0,0,0,0,0,0,115,117,112,101,0,0,0,0,99,108,111,115,101,112,97,116,104,32,115,116,114,111,107,101,10,0,0,0,0,0,0,0,47,103,110,98,117,53,47,50,0,0,0,0,0,0,0,0,72,101,108,118,101,116,105,99,97,0,0,0,0,0,0,0,37,100,44,37,100,44,37,100,44,37,100,0,0,0,0,0,35,32,71,101,110,101,114,97,116,101,100,32,98,121,32,37,115,32,118,101,114,115,105,111,110,32,37,115,32,40,37,115,41,10,0,0,0,0,0,0,49,46,50,0,0,0,0,0,112,110,103,58,100,111,116,0,33,40,113,45,62,113,116,115,41,0,0,0,0,0,0,0,115,117,112,51,0,0,0,0,46,46,47,46,46,47,46,46,47,108,105,98,47,110,101,97,116,111,103,101,110,47,110,101,97,116,111,105,110,105,116,46,99,0,0,0,0,0,0,0,99,111,109,112,114,101,115,115,0,0,0,0,0,0,0,0,47,103,110,98,117,53,47,49,0,0,0,0,0,0,0,0,97,103,116,97,105,108,40,101,41,32,61,61,32,85,70,95,102,105,110,100,40,97,103,116,97,105,108,40,101,41,41,0,47,97,99,99,101,110,116,51,47,49,0,0,0,0,0,0,115,117,112,50,0,0,0,0,97,108,115,101,0,0,0,0,104,101,120,97,103,111,110,0,37,46,53,103,0,0,0,0,115,116,114,105,99,116,0,0,47,103,110,98,117,52,47,52,0,0,0,0,0,0,0,0,116,0,0,0,114,0,0,0,117,0,0,0,101,0,0,0,0,0,0,0,0,0,0,0,115,117,112,49,0,0,0,0,47,103,110,98,117,52,47,51,0,0,0,0,0,0,0,0,115,117,112,0,0,0,0,0,116,101,101,0,0,0,0,0,47,103,110,98,117,52,47,50,0,0,0,0,0,0,0,0,115,117,109,0,0,0,0,0,47,103,110,98,117,52,47,49,0,0,0,0,0,0,0,0,97,103,108,115,116,110,111,100,101,0,0,0,0,0,0,0,115,117,98,101,0,0,0,0,99,104,97,114,0,0,0,0,47,103,110,98,117,51,47,51,0,0,0,0,0,0,0,0,115,117,98,0,0,0,0,0,47,103,110,98,117,51,47,50,0,0,0,0,0,0,0,0,47,103,110,98,117,51,47,49,0,0,0,0,0,0,0,0,38,103,116,59,0,0,0,0,115,112,97,100,101,115,0,0,35,32,0,0,0,0,0,0,87,97,114,110,105,110,103,0,47,97,99,99,101,110,116,55,47,54,0,0,0,0,0,0,115,105,109,0,0,0,0,0,47,100,97,114,107,50,56,47,56,0,0,0,0,0,0,0,121,101,108,108,111,119,0,0,32,99,114,101,97,116,101,32,111,118,97,108,32,0,0,0,87,69,66,80,0,0,0,0,100,97,114,107,103,111,108,100,101,110,114,111,100,0,0,0,115,105,103,109,97,102,0,0,99,108,111,115,101,112,97,116,104,32,102,105,108,108,10,0,100,97,114,107,112,117,114,112,108,101,0,0,0,0,0,0,47,100,97,114,107,50,56,47,55,0,0,0,0,0,0,0,67,111,117,114,105,101,114,45,66,111,108,100,79,98,108,105,113,117,101,0,0,0,0,0,37,100,44,37,100,44,37,100,0,0,0,0,0,0,0,0,35,70,73,71,32,51,46,50,10,0,0,0,0,0,0,0,120,100,111,116,118,101,114,115,105,111,110,0,0,0,0,0,115,118,103,58,109,97,112,0,97,108,105,99,101,98,108,117,101,0,0,0,0,0,0,0,115,101,110,100,32,114,97,110,100,111,109,32,99,111,111,114,100,105,110,97,116,101,115,10,0,0,0,0,0,0,0,0,113,45,62,110,32,61,61,32,49,0,0,0,0,0,0,0,115,105,103,109,97,0,0,0,115,99,97,108,105,110,103,0,47,100,97,114,107,50,56,47,54,0,0,0,0,0,0,0,97,103,104,101,97,100,40,101,41,32,61,61,32,85,70,95,102,105,110,100,40,97,103,104,101,97,100,40,101,41,41,0,115,104,121,0,0,0,0,0,111,109,112,111,117,110,100,0,112,101,110,116,97,103,111,110,0,0,0,0,0,0,0,0,75,0,0,0,0,0,0,0,44,37,46,53,103,0,0,0,93,59,10,0,0,0,0,0,47,100,97,114,107,50,56,47,53,0,0,0,0,0,0,0,112,101,110,99,111,108,111,114,0,0,0,0,0,0,0,0,112,115,50,58,112,115,0,0,115,101,99,116,0,0,0,0,47,100,97,114,107,50,56,47,52,0,0,0,0,0,0,0,115,100,111,116,0,0,0,0,47,100,97,114,107,50,56,47,51,0,0,0,0,0,0,0,112,111,118,0,0,0,0,0,115,99,97,114,111,110,0,0,47,100,97,114,107,50,56,47,50,0,0,0,0,0,0,0,97,103,112,114,118,110,111,100,101,0,0,0,0,0,0,0,47,37,115,47,37,115,0,0,37,115,32,37,115,10,0,0,88,49,49,47,0,0,0,0,105,103,104,116,103,114,101,121,0,0,0,0,0,0,0,0,104,105,116,101,0,0,0,0,115,98,113,117,111,0,0,0,108,97,99,107,0,0,0,0,99,109,97,112,58,109,97,112,0,0,0,0,0,0,0,0,121,101,108,108,111,119,52,0,121,101,108,108,111,119,51,0,47,100,97,114,107,50,56,47,49,0,0,0,0,0,0,0,121,101,108,108,111,119,50,0,121,101,108,108,111,119,49,0,119,104,101,97,116,52,0,0,114,115,113,117,111,0,0,0,119,104,101,97,116,51,0,0,119,104,101,97,116,50,0,0,119,104,101,97,116,49,0,0,47,100,97,114,107,50,55,47,55,0,0,0,0,0,0,0,118,105,111,108,101,116,114,101,100,52,0,0,0,0,0,0,118,105,111,108,101,116,114,101,100,51,0,0,0,0,0,0,118,105,111,108,101,116,114,101,100,50,0,0,0,0,0,0,118,105,111,108,101,116,114,101,100,49,0,0,0,0,0,0,106,112,101,103,58,115,118,103,0,0,0,0,0,0,0,0,115,102,100,112,0,0,0,0,114,115,97,113,117,111,0,0,116,117,114,113,117,111,105,115,101,52,0,0,0,0,0,0,116,117,114,113,117,111,105,115,101,51,0,0,0,0,0,0,116,117,114,113,117,111,105,115,101,50,0,0,0,0,0,0,98,111,120,0,0,0,0,0,116,117,114,113,117,111,105,115,101,49,0,0,0,0,0,0,47,100,97,114,107,50,55,47,54,0,0,0,0,0,0,0,85,115,105,110,103,32,100,101,102,97,117,108,116,32,99,97,108,99,117,108,97,116,105,111,110,32,102,111,114,32,114,111,111,116,32,110,111,100,101,10,0,0,0,0,0,0,0,0,83,112,97,114,115,101,77,97,116,114,105,120,95,105,115,95,115,121,109,109,101,116,114,105,99,40,66,44,32,70,65,76,83,69,41,0,0,0,0,0,116,111,109,97,116,111,52,0,116,111,109,97,116,111,51,0,47,97,99,99,101,110,116,55,47,53,0,0,0,0,0,0,116,111,109,97,116,111,50,0,116,111,109,97,116,111,49,0,114,108,109,0,0,0,0,0,37,115,58,32,0,0,0,0,116,104,105,115,116,108,101,52,0,0,0,0,0,0,0,0,116,104,105,115,116,108,101,51,0,0,0,0,0,0,0,0,116,104,105,115,116,108,101,50,0,0,0,0,0,0,0,0,47,100,97,114,107,50,55,47,53,0,0,0,0,0,0,0,108,101,110,32,62,32,48,0,116,104,105,115,116,108,101,49,0,0,0,0,0,0,0,0,119,104,105,116,101,0,0,0,32,45,111,117,116,108,105,110,101,32,0,0,0,0,0,0,116,97,110,52,0,0,0,0,116,97,110,51,0,0,0,0,60,33,45,45,32,71,101,110,101,114,97,116,101,100,32,98,121,32,0,0,0,0,0,0,116,97,110,50,0,0,0,0,100,97,114,107,99,121,97,110,0,0,0,0,0,0,0,0,103,97,105,110,32,60,61,32,113,45,62,110,103,97,105,110,0,0,0,0,0,0,0,0,116,97,110,49,0,0,0,0,114,104,111,0,0,0,0,0,32,99,117,114,118,101,116,111,10,0,0,0,0,0,0,0,115,116,101,101,108,98,108,117,101,52,0,0,0,0,0,0,115,116,101,101,108,98,108,117,101,51,0,0,0,0,0,0,115,116,101,101,108,98,108,117,101,50,0,0,0,0,0,0,47,100,97,114,107,50,55,47,52,0,0,0,0,0,0,0,67,111,117,114,105,101,114,0,101,112,115,105,108,111,110,0,32,99,111,111,114,100,115,61,34,0,0,0,0,0,0,0,35,32,101,110,100,32,111,102,32,70,73,71,32,102,105,108,101,10,0,0,0,0,0,0,115,116,101,101,108,98,108,117,101,49,0,0,0,0,0,0,95,100,114,97,119,95,0,0,101,112,115,58,109,97,112,0,103,114,97,112,104,32,105,115,32,100,105,115,99,111,110,110,101,99,116,101,100,46,32,72,101,110,99,101,44,32,116,104])
.concat([101,32,99,105,114,99,117,105,116,32,109,111,100,101,108,10,0,0,0,0,0,0,0,0,115,112,114,105,110,103,103,114,101,101,110,52,0,0,0,0,100,109,101,97,110,32,61,32,37,102,44,32,114,104,111,32,61,32,37,102,10,0,0,0,115,112,114,105,110,103,103,114,101,101,110,51,0,0,0,0,113,45,62,113,116,115,91,105,105,93,0,0,0,0,0,0,32,34,104,116,116,112,58,47,47,119,119,119,46,119,51,46,111,114,103,47,71,114,97,112,104,105,99,115,47,83,86,71,47,49,46,49,47,68,84,68,47,115,118,103,49,49,46,100,116,100,34,62,10,0,0,0,115,112,114,105,110,103,103,114,101,101,110,50,0,0,0,0,115,112,114,105,110,103,103,114,101,101,110,49,0,0,0,0,114,102,108,111,111,114,0,0,102,97,105,108,101,100,32,116,111,32,105,110,105,116,32,108,105,98,108,116,100,108,10,0,109,97,107,101,65,100,100,80,111,108,121,58,32,117,110,107,110,111,119,110,32,115,104,97,112,101,32,116,121,112,101,32,37,115,10,0,0,0,0,0,115,110,111,119,52,0,0,0,115,110,111,119,51,0,0,0,115,110,111,119,50,0,0,0,47,100,97,114,107,50,55,47,51,0,0,0,0,0,0,0,40,78,68,95,85,70,95,115,105,122,101,40,110,41,32,60,61,32,49,41,32,124,124,32,40,110,32,61,61,32,108,101,97,100,101,114,41,0,0,0,102,108,97,116,105,110,100,101,120,40,97,103,116,97,105,108,40,101,41,41,32,60,32,77,45,62,110,99,111,108,115,0,115,110,111,119,49,0,0,0,46,46,47,46,46,47,46,46,47,108,105,98,47,110,101,97,116,111,103,101,110,47,110,101,97,116,111,115,112,108,105,110,101,115,46,99,0,0,0,0,115,108,97,116,101,103,114,97,121,52,0,0,0,0,0,0,60,33,68,79,67,84,89,80,69,32,115,118,103,32,80,85,66,76,73,67,32,34,45,47,47,87,51,67,47,47,68,84,68,32,83,86,71,32,49,46,49,47,47,69,78,34,10,0,115,108,97,116,101,103,114,97,121,51,0,0,0,0,0,0,115,108,97,116,101,103,114,97,121,50,0,0,0,0,0,0,114,101,103,0,0,0,0,0,117,114,118,101,100,0,0,0,115,108,97,116,101,103,114,97,121,49,0,0,0,0,0,0,76,97,121,111,117,116,32,116,121,112,101,58,32,34,37,115,34,32,110,111,116,32,114,101,99,111,103,110,105,122,101,100,46,32,85,115,101,32,111,110,101,32,111,102,58,37,115,10,0,0,0,0,0,0,0,0,104,111,117,115,101,0,0,0,37,108,102,44,37,108,102,37,99,0,0,0,0,0,0,0,115,108,97,116,101,98,108,117,101,52,0,0,0,0,0,0,37,46,53,103,44,37,46,53,103,44,37,46,53,103,0,0,115,108,97,116,101,98,108,117,101,51,0,0,0,0,0,0,47,100,97,114,107,50,55,47,50,0,0,0,0,0,0,0,115,108,97,116,101,98,108,117,101,50,0,0,0,0,0,0,115,108,97,116,101,98,108,117,101,49,0,0,0,0,0,0,108,105,103,104,116,0,0,0,115,107,121,98,108,117,101,52,0,0,0,0,0,0,0,0,34,32,116,121,112,101,61,34,116,101,120,116,47,99,115,115,34,63,62,10,0,0,0,0,115,107,121,98,108,117,101,51,0,0,0,0,0,0,0,0,115,107,121,98,108,117,101,50,0,0,0,0,0,0,0,0,114,101,97,108,0,0,0,0,115,107,121,98,108,117,101,49,0,0,0,0,0,0,0,0,115,105,101,110,110,97,52,0,115,105,101,110,110,97,51,0,47,100,97,114,107,50,55,47,49,0,0,0,0,0,0,0,115,105,101,110,110,97,50,0,115,105,101,110,110,97,49,0,115,101,97,115,104,101,108,108,52,0,0,0,0,0,0,0,100,101,108,97,117,110,97,121,95,116,114,105,58,32,37,115,10,0,0,0,0,0,0,0,60,63,120,109,108,45,115,116,121,108,101,115,104,101,101,116,32,104,114,101,102,61,34,0,115,101,97,115,104,101,108,108,51,0,0,0,0,0,0,0,98,101,115,116,99,111,115,116,32,60,32,72,85,71,69,95,86,65,76,0,0,0,0,0,115,101,97,115,104,101,108,108,50,0,0,0,0,0,0,0,114,100,113,117,111,0,0,0,115,101,97,115,104,101,108,108,49,0,0,0,0,0,0,0,115,101,97,103,114,101,101,110,52,0,0,0,0,0,0,0,115,101,97,103,114,101,101,110,51,0,0,0,0,0,0,0,47,100,97,114,107,50,54,47,54,0,0,0,0,0,0,0,115,101,97,103,114,101,101,110,50,0,0,0,0,0,0,0,115,101,97,103,114,101,101,110,49,0,0,0,0,0,0,0,115,116,121,108,101,115,104,101,101,116,0,0,0,0,0,0,84,48,0,0,0,0,0,0,115,97,108,109,111,110,52,0,115,97,108,109,111,110,51,0,114,99,101,105,108,0,0,0,115,97,108,109,111,110,50,0,104,101,105,103,104,116,0,0,115,97,108,109,111,110,49,0,99,99,37,115,43,37,100,0,116,111,111,32,109,97,110,121,32,40,62,32,37,100,41,32,115,97,109,101,123,104,101,97,100,44,116,97,105,108,125,32,103,114,111,117,112,115,32,102,111,114,32,110,111,100,101,32,37,115,10,0,0,0,0,0,47,100,97,114,107,50,54,47,53,0,0,0,0,0,0,0,97,103,110,120,116,110,111,100,101,0,0,0,0,0,0,0,114,111,121,97,108,98,108,117,101,52,0,0,0,0,0,0,114,111,121,97,108,98,108,117,101,51,0,0,0,0,0,0,114,111,121,97,108,98,108,117,101,50,0,0,0,0,0,0,114,111,121,97,108,98,108,117,101,49,0,0,0,0,0,0,60,63,120,109,108,32,118,101,114,115,105,111,110,61,34,49,46,48,34,32,101,110,99,111,100,105,110,103,61,34,85,84,70,45,56,34,32,115,116,97,110,100,97,108,111,110,101,61,34,110,111,34,63,62,10,0,114,111,115,121,98,114,111,119,110,52,0,0,0,0,0,0,75,80,95,76,101,102,116,0,114,97,114,114,0,0,0,0,114,111,115,121,98,114,111,119,110,51,0,0,0,0,0,0,114,111,115,121,98,114,111,119,110,50,0,0,0,0,0,0,114,111,115,121,98,114,111,119,110,49,0,0,0,0,0,0,46,46,47,46,46,47,46,46,47,108,105,98,47,100,111,116,103,101,110,47,109,105,110,99,114,111,115,115,46,99,0,0,71,68,95,109,105,110,114,97,110,107,40,103,41,32,61,61,32,48,0,0,0,0,0,0,47,100,97,114,107,50,54,47,52,0,0,0,0,0,0,0,114,101,100,52,0,0,0,0,114,101,100,51,0,0,0,0,110,32,33,61,32,78,68,95,110,101,120,116,40,110,41,0,114,101,100,50,0,0,0,0,114,101,100,49,0,0,0,0,32,120,109,108,110,115,58,120,108,105,110,107,61,34,104,116,116,112,58,47,47,119,119,119,46,119,51,46,111,114,103,47,49,57,57,57,47,120,108,105,110,107,34,0,0,0,0,0,69,68,95,108,97,98,101,108,40,102,101,41,0,0,0,0,112,117,114,112,108,101,52,0,114,97,113,117,111,0,0,0,112,117,114,112,108,101,51,0,99,111,109,112,111,117,110,100,0,0,0,0,0,0,0,0,100,101,103,101,110,101,114,97,116,101,32,99,111,110,99,101,110,116,114,97,116,101,100,32,114,97,110,107,32,37,115,44,37,100,10,0,0,0,0,0,112,117,114,112,108,101,50,0,112,117,114,112,108,101,49,0,37,115,32,45,62,32,37,115,58,32,115,112,108,105,110,101,32,115,105,122,101,32,62,32,49,32,110,111,116,32,115,117,112,112,111,114,116,101,100,10,0,0,0,0,0,0,0,0,78,68,95,114,97,110,107,40,102,114,111,109,41,32,60,32,78,68,95,114,97,110,107,40,116,111,41,0,0,0,0,0,47,100,97,114,107,50,54,47,51,0,0,0,0,0,0,0,69,68,95,116,111,95,118,105,114,116,40,111,114,105,103,41,32,61,61,32,78,85,76,76,0,0,0,0,0,0,0,0,112,108,117,109,52,0,0,0,78,111,32,108,105,98,122,32,115,117,112,112,111,114,116,46,10,0,0,0,0,0,0,0,112,108,117,109,51,0,0,0,112,108,117,109,50,0,0,0,32,120,109,108,110,115,61,34,104,116,116,112,58,47,47,119,119,119,46,119,51,46,111,114,103,47,50,48,48,48,47,115,118,103,34,0,0,0,0,0,112,108,117,109,49,0,0,0,109,105,110,100,105,115,116,0,114,97,110,103,0,0,0,0,112,105,110,107,52,0,0,0,112,105,110,107,51,0,0,0,112,105,110,107,50,0,0,0,50,48,49,51,48,54,49,54,46,48,48,49,49,0,0,0,112,105,110,107,49,0,0,0,47,100,97,114,107,50,54,47,50,0,0,0,0,0,0,0,115,111,117,114,99,101,32,112,111,105,110,116,32,110,111,116,32,105,110,32,97,110,121,32,116,114,105,97,110,103,108,101,0,0,0,0,0,0,0,0,112,101,97,99,104,112,117,102,102,52,0,0,0,0,0,0,112,101,97,99,104,112,117,102,102,51,0,0,0,0,0,0,32,118,105,101,119,66,111,120,61,34,37,46,50,102,32,37,46,50,102,32,37,46,50,102,32,37,46,50,102,34,0,0,47,97,99,99,101,110,116,55,47,52,0,0,0,0,0,0,112,101,97,99,104,112,117,102,102,50,0,0,0,0,0,0,112,101,97,99,104,112,117,102,102,49,0,0,0,0,0,0,114,97,100,105,99,0,0,0,112,97,108,101,118,105,111,108,101,116,114,101,100,52,0,0,112,97,108,101,118,105,111,108,101,116,114,101,100,51,0,0,47,100,97,114,107,50,54,47,49,0,0,0,0,0,0,0,37,102,0,0,0,0,0,0,112,97,108,101,118,105,111,108,101,116,114,101,100,50,0,0,112,97,108,101,118,105,111,108,101,116,114,101,100,49,0,0,116,101,97,108,0,0,0,0,112,97,108,101,116,117,114,113,117,111,105,115,101,52,0,0,60,115,118,103,32,119,105,100,116,104,61,34,37,100,112,116,34,32,104,101,105,103,104,116,61,34,37,100,112,116,34,10,0,0,0,0,0,0,0,0,60,115,118,103,0,0,0,0,112,97,108,101,116,117,114,113,117,111,105,115,101,51,0,0,100,97,114,107,98,108,117,101,0,0,0,0,0,0,0,0,112,97,108,101,116,117,114,113,117,111,105,115,101,50,0,0,114,65,114,114,0,0,0,0,112,97,108,101,116,117,114,113,117,111,105,115,101,49,0,0,115,116,114,111,107,101,10,0,112,97,108,101,103,114,101,101,110,52,0,0,0,0,0,0,112,97,108,101,103,114,101,101,110,51,0,0,0,0,0,0,47,100,97,114,107,50,53,47,53,0,0,0,0,0,0,0,67,111,117,114,105,101,114,45,66,111,108,100,0,0,0,0,32,97,108,116,61,34,34,0,112,97,108,101,103,114,101,101,110,50,0,0,0,0,0,0,112,97,108,101,103,114,101,101,110,49,0,0,0,0,0,0,112,115,58,109,97,112,0,0,105,105,32,60,32,49,60,60,100,105,109,32,38,38,32,105,105,32,62,61,32,48,0,0,32,80,97,103,101,115,58,32,37,100,32,45,45,62,10,0,98,97,100,32,101,100,103,101,32,108,101,110,32,34,37,115,34,0,0,0,0,0,0,0,111,114,99,104,105,100,52,0,111,114,99,104,105,100,51,0,113,117,111,116,0,0,0,0,116,104,101,32,103,114,97,112,104,32,105,110,116,111,32,99,111,110,110,101,99,116,101,100,32,99,111,109,112,111,110,101,110,116,115,46,10,0,0,0,111,114,99,104,105,100,50,0,111,114,99,104,105,100,49,0,86,111,114,111,110,111,105,0,103,114,97,112,104,32,37,115,44,32,99,111,111,114,100,32,37,115,44,32,101,120,112,101,99,116,101,100,32,102,111,117,114,32,100,111,117,98,108,101,115,10,0,0,0,0,0,0,111,114,97,110,103,101,114,101,100,52,0,0,0,0,0,0,47,100,97,114,107,50,53,47,52,0,0,0,0,0,0,0,108,101,97,100,101,114,32,33,61,32,78,85,76,76,0,0,102,108,97,116,105,110,100,101,120,40,97,103,104,101,97,100,40,101,41,41,32,60,32,77,45,62,110,114,111,119,115,0,111,114,97,110,103,101,114,101,100,51,0,0,0,0,0,0,111,114,97,110,103,101,114,101,100,50,0,0,0,0,0,0,111,114,97,110,103,101,114,101,100,49,0,0,0,0,0,0,111,114,97,110,103,101,52,0,111,114,97,110,103,101,51,0,112,115,105,0,0,0,0,0,85,84,70,56,32,99,111,100,101,115,32,62,32,51,32,98,121,116,101,115,32,97,114,101,32,110,111,116,32,99,117,114,114,101,110,116,108,121,32,115,117,112,112,111,114,116,101,100,32,40,103,114,97,112,104,32,37,115,41,32,45,32,116,114,101,97,116,101,100,32,97,115,32,76,97,116,105,110,45,49,46,32,80,101,114,104,97,112,115,32,34,45,71,99,104,97,114,115,101,116,61,108,97,116,105,110,49,34,32,105,115,32,110,101,101,100,101,100,63,10,0,0,0,0,0,0,0,0,111,114,97,110,103,101,50,0,111,114,97,110,103,101,49,0,112,97,114,97,108,108,101,108,111,103,114,97,109,0,0,0,111,108,105,118,101,100,114,97,98,52,0,0,0,0,0,0,117,32,61,61,32,85,70,95,102,105,110,100,40,117,41,0,47,112,114,111,99,47,115,101,108,102,47,109,97,112,115,0,47,100,97,114,107,50,53,47,51,0,0,0,0,0,0,0,111,108,105,118,101,100,114,97,98,51,0,0,0,0,0,0,80,45,62,101,110,100,46,116,104,101,116,97,32,60,32,50,32,42,32,77,95,80,73,0,111,108,105,118,101,100,114,97,98,50,0,0,0,0,0,0,111,108,105,118,101,100,114,97,98,49,0,0,0,0,0,0,60,33,45,45,0,0,0,0,112,114,111,112,0,0,0,0,110,97,118,97,106,111,119,104,105,116,101,52,0,0,0,0,110,97,118,97,106,111,119,104,105,116,101,51,0,0,0,0,47,100,97,114,107,50,53,47,50,0,0,0,0,0,0,0,110,97,118,97,106,111,119,104,105,116,101,50,0,0,0,0,110,97,118,97,106,111,119,104,105,116,101,49,0,0,0,0,115,105,100,101,115,32,61,61,32,52,0,0,0,0,0,0,60,47,115,118,103,62,10,0,109,105,115,116,121,114,111,115,101,52,0,0,0,0,0,0,99,97,110,110,111,116,32,114,101,45,97,108,108,111,99,97,116,101,32,112,115,10,0,0,109,105,115,116,121,114,111,115,101,51,0,0,0,0,0,0,112,114,111,100,0,0,0,0,109,105,115,116,121,114,111,115,101,50,0,0,0,0,0,0,109,105,115,116,121,114,111,115,101,49,0,0,0,0,0,0,108,111,115,116,32,37,115,32,37,115,32,101,100,103,101,10,0,0,0,0,0,0,0,0,47,100,97,114,107,50,53,47,49,0,0,0,0,0,0,0,110,111,100,101,32,0,0,0,34,32,99,108,97,115,115,61,34,108,97,121,101,114,34,62,10,0,0,0,0,0,0,0,112,114,105,109,101,0,0,0,109,101,100,105,117,109,112,117,114,112,108,101,52,0,0,0,109,101,100,105,117,109,112,117,114,112,108,101,51,0,0,0,105,110,32,108,97,98,101,108,32,111,102,32,101,100,103,101,32,37,115,32,37,115,32,37,115,10,0,0,0,0,0,0,109,101,100,105,117,109,112,117,114,112,108,101,50,0,0,0,109,101,100,105,117,109,112,117,114,112,108,101,49,0,0,0,47,100,97,114,107,50,52,47,52,0,0,0,0,0,0,0,97,103,102,115,116,110,111,100,101,0,0,0,0,0,0,0,109,101,100,105,117,109,111,114,99,104,105,100,52,0,0,0,109,101,100,105,117,109,111,114,99,104,105,100,51,0,0,0,109,101,100,105,117,109,111,114,99,104,105,100,50,0,0,0,32,116,114,97,110,115,102,111,114,109,61,34,115,99,97,108,101,40,37,103,32,37,103,41,32,114,111,116,97,116,101,40,37,100,41,32,116,114,97,110,115,108,97,116,101,40,37,103,32,37,103,41,34,62,10,0,109,101,100,105,117,109,111,114,99,104,105,100,49,0,0,0,112,111,117,110,100,0,0,0,109,97,114,111,111,110,52,0,109,97,114,111,111,110,51,0,47,100,97,114,107,50,52,47,51,0,0,0,0,0,0,0,109,97,114,111,111,110,50,0,109,97,114,111,111,110,49,0,109,97,103,101,110,116,97,52,0,0,0,0,0,0,0,0,34,32,99,108,97,115,115,61,34,103,114,97,112,104,34,0,109,97,103,101,110,116,97,51,0,0,0,0,0,0,0,0,109,97,103,101,110,116,97,50,0,0,0,0,0,0,0,0,112,108,117,115,109,110,0,0,109,97,103,101,110,116,97,49,0,0,0,0,0,0,0,0,47,100,97,114,107,50,52,47,50,0,0,0,0,0,0,0,115,121,110,116,97,120,32,101,114,114,111,114,0,0,0,0,108,105,103,104,116,121,101,108,108,111,119,52,0,0,0,0,108,105,103,104,116,121,101,108,108,111,119,51,0,0,0,0,108,105,103,104,116,121,101,108,108,111,119,50,0,0,0,0,108,105,103,104,116,121,101,108,108,111,119,49,0,0,0,0,34,32,99,108,97,115,115,61,34,99,108,117,115,116,101,114,34,62,0,0,0,0,0,0,108,105,103,104,116,115,116,101,101,108,98,108,117,101,52,0,112,105,118,0,0,0,0,0,108,105,103,104,116,115,116,101,101,108,98,108,117,101,51,0,108,105,103,104,116,115,116,101,101,108,98,108,117,101,50,0,108,105,103,104,116,115,116,101,101,108,98,108,117,101,49,0,47,100,97,114,107,50,52,47,49,0,0,0,0,0,0,0,108,105,103,104,116,115,108,97,116,101,98,108,117,101,0,0,108,105,103,104,116,115,107,121,98,108,117,101,52,0,0,0,34,32,99,108,97,115,115,61,34,110,111,100,101,34,62,0,47,97,99,99,101,110,116,55,47,51,0,0,0,0,0,0,108,105,103,104,116,115,107,121,98,108,117,101,51,0,0,0,108,105,103,104,116,115,107,121,98,108,117,101,50,0,0,0,112,105,0,0,0,0,0,0,99,111,117,114,0,0,0,0,108,105,103,104,116,115,107,121,98,108,117,101,49,0,0,0,108,105,103,104,116,115,97,108,109,111,110,52,0,0,0,0,47,100,97,114,107,50,51,47,51,0,0,0,0,0,0,0,108,105,103,104,116,115,97,108,109,111,110,51,0,0,0,0,108,105,103,104,116,115,97,108,109,111,110,50,0,0,0,0,115,105,108,118,101,114,0,0,32,99,114,101,97,116,101,32,112,111,108,121,103,111,110,32,0,0,0,0,0,0,0,0,108,105,103,104,116,115,97,108,109,111,110,49,0,0,0,0,95,37,115,0,0,0,0,0,108,105,103,104,116,112,105,110,107,52,0,0,0,0,0,0,99,121,97,110,0,0,0,0,108,105,103,104,116,112,105,110,107,51,0,0,0,0,0,0,112,104,105,0,0,0,0,0,108,105,103,104,116,112,105,110,107,50,0,0,0,0,0,0,32,108,105,110,101,116,111,10,0,0,0,0,0,0,0,0,108,105,103,104,116,112,105,110,107,49,0,0,0,0,0,0,32,32,34,37,115,34,10,0,47,100,97,114,107,50,51,47,50,0,0,0,0,0,0,0,84,105,109,101,115,45,66,111,108,100,73,116,97,108,105,99,0,0,0,0,0,0,0,0,32,116,105,116,108,101,61,34,0,0,0,0,0,0,0,0,37,100,32,37,100,32,37,100,32,37,100,32,37,100,32,37,100,32,37,46,49,102,32,37,46,52,102,32,37,100,32,37,46,49,102,32,37,46,49,102,32,37,100,32,37,100,32,37,115,92,48,48,49,10,0,0,70,32,37,102,32,0,0,0,106,112,103,58,109,97,112,0,108,105,103,104,116,103,111,108,100,101,110,114,111,100,52,0,115,112,114,105,110,103,95,101,108,101,99,116,114,105,99,97,108,95,101,109,98,101,100,100,105,110,103,95,115,108,111,119,0,0,0,0,0,0,0,0,100,105,109,0,0,0,0,0,37,112,0,0,0,0,0,0,108,105,103,104,116,103,111,108,100,101,110,114,111,100,51,0,33,40,113,45,62,108,41,0,60,47,116,105,116,108,101,62,10,0,0,0,0,0,0,0,37,108,102,0,0,0,0,0,108,105,103,104,116,103,111,108,100,101,110,114,111,100,50,0,108,105,103,104,116,103,111,108,100,101,110,114,111,100,49,0,112,101,114,112,0,0,0,0,65,108,116,101,114,110,97,116,105,118,101,108,121,44,32,99,111,110,115,105,100,101,114,32,114,117,110,110,105,110,103,32,110,101,97,116,111,32,117,115,105,110,103,32,45,71,112,97,99,107,61,116,114,117,101,32,111,114,32,100,101,99,111,109,112,111,115,105,110,103,10,0,108,105,103,104,116,103,111,108,100,101,110,114,111,100,0,0,108,105,103,104,116,99,121,97,110,52,0,0,0,0,0,0,107,101,121,0,0,0,0,0,108,105,103,104,116,99,121,97,110,51,0,0,0,0,0,0,118,111,114,111,110,111,105,0,37,108,102,44,37,108,102,44,37,108,102,44,37,108,102,37,99,0,0,0,0,0,0,0,108,105,103,104,116,99,121,97,110,50,0,0,0,0,0,0,47,100,97,114,107,50,51,47,49,0,0,0,0,0,0,0,46,46,47,46,46,47,46,46,47,108,105,98,47,100,111,116,103,101,110,47,114,97,110,107,46,99,0,0,0,0,0,0,78,68,95,114,97,110,107,40,118,41,32,61,61,32,114,0,108,105,103,104,116,99,121,97,110,49,0,0,0,0,0,0,47,98,117,112,117,57,47,57,0,0,0,0,0,0,0,0,108,105,103,104,116,98,108,117,101,52,0,0,0,0,0,0,92,69,0,0,0,0,0,0,108,105,103,104,116,98,108,117,101,51,0,0,0,0,0,0,108,105,103,104,116,98,108,117,101,50,0,0,0,0,0,0,112,101,114,109,105,108,0,0,73,110,118,97,108,105,100,32,51,45,98,121,116,101,32,85,84,70,56,32,102,111,117,110,100,32,105,110,32,105,110,112,117,116,32,111,102,32,103,114,97,112,104,32,37,115,32,45,32,116,114,101,97,116,101,100,32,97,115,32,76,97,116,105,110,45,49,46,32,80,101,114,104,97,112,115,32,34,45,71,99,104,97,114,115,101,116,61,108,97,116,105,110,49,34,32,105,115,32,110,101,101,100,101,100,63,10,0,0,0,0,0,108,105,103,104,116,98,108,117,101,49,0,0,0,0,0,0,116,114,97,112,101,122,105,117,109,0,0,0,0,0,0,0,108,101,109,111,110,99,104,105,102,102,111,110,52,0,0,0,37,46,50,102,0,0,0,0,108,101,109,111,110,99,104,105,102,102,111,110,51,0,0,0,108,101,109,111,110,99,104,105,102,102,111,110,50,0,0,0,108,101,109,111,110,99,104,105,102,102,111,110,49,0,0,0,108,97,118,101,110,100,101,114,98,108,117,115,104,52,0,0,60,116,105,116,108,101,62,0,102,97,116,97,108,32,102,108,101,120,32,115,99,97,110,110,101,114,32,105,110,116,101,114,110,97,108,32,101,114,114,111,114,45,45,110,111,32,97,99,116,105,111,110,32,102,111,117,110,100,0,0,0,0,0,0,108,97,118,101,110,100,101,114,98,108,117,115,104,51,0,0,112,97,114,116,0,0,0,0,108,97,118,101,110,100,101,114,98,108,117,115,104,50,0,0,108,97,118,101,110,100,101,114,98,108,117,115,104,49,0,0,47,98,117,112,117,57,47,56,0,0,0,0,0,0,0,0,107,104,97,107,105,52,0,0,107,104,97,107,105,51,0,0,107,104,97,107,105,50,0,0,37,115,10,0,0,0,0,0,107,104,97,107,105,49,0,0,34,32,99,108,97,115,115,61,34,101,100,103,101,34,62,0,105,118,111,114,121,52,0,0,112,97,114,97,0,0,0,0,105,118,111,114,121,51,0,0,105,118,111,114,121,50,0,0,105,118,111,114,121,49,0,0,47,98,117,112,117,57,47,55,0,0,0,0,0,0,0,0,105,110,100,105,97,110,114,101,100,52,0,0,0,0,0,0,105,110,100,105,97,110,114,101,100,51,0,0,0,0,0,0,60,103,32,105,100,61,34,0,105,110,100,105,97,110,114,101,100,50,0,0,0,0,0,0,105,110,100,105,97,110,114,101,100,49,0,0,0,0,0,0,111,117,109,108,0,0,0,0,104,111,116,112,105,110,107,52,0,0,0,0,0,0,0,0,104,111,116,112,105,110,107,51,0,0,0,0,0,0,0,0,104,111,116,112,105,110,107,50,0,0,0,0,0,0,0,0,47,98,117,112,117,57,47,54,0,0,0,0,0,0,0,0,97,103,105,100,110,111,100,101,0,0,0,0,0,0,0,0,104,111,116,112,105,110,107,49,0,0,0,0,0,0,0,0,104,111,110,101,121,100,101,119,52,0,0,0,0,0,0,0,104,111,110,101,121,100,101,119,51,0,0,0,0,0,0,0,104,111,110,101,121,100,101,119,50,0,0,0,0,0,0,0,104,111,110,101,121,100,101,119,49,0,0,0,0,0,0,0,111,116,105,109,101,115,0,0,103,114,101,121,57,57,0,0,103,114,101,121,57,56,0,0,103,114,101,121,57,55,0,0,47,98,117,112,117,57,47,53,0,0,0,0,0,0,0,0,103,114,101,121,57,54,0,0,103,114,101,121,57,53,0,0,103,114,101,121,57,52,0,0,103,114,101,121,57,51,0,0,32,116,97,114,103,101,116,61,34,0,0,0,0,0,0,0,103,114,101,121,57,50,0,0,103,114,101,121,57,49,0,0,111,116,105,108,100,101,0,0,103,114,101,121,57,48,0,0,103,114,101,121,57,0,0,0,103,114,101,121,56,57,0,0,103,114,101,121,56,56,0,0,97,103,114,111,111,116,32,111,102,32,97,32,98,97,100,32,111,98,106,101,99,116,0,0,47,98,117,112,117,57,47,52,0,0,0,0,0,0,0,0,103,114,101,121,56,55,0,0,103,114,101,121,56,54,0,0,103,114,101,121,56,53,0,0,103,114,101,121,56,52,0,0,32,120,108,105,110,107,58,116,105,116,108,101,61,34,0,0,103,114,101,121,56,51,0,0,103,114,101,121,56,50,0,0,111,115,108,97,115,104,0,0,103,114,101,121,56,49,0,0,103,114,101,121,56,48,0,0,103,114,101,121,56,0,0,0,103,114,101,121,55,57,0,0,47,98,117,112,117,57,47,51,0,0,0,0,0,0,0,0,103,114,101,121,55,56,0,0,103,114,101,121,55,55,0,0,103,114,101,121,55,54,0,0,103,114,101,121,55,53,0,0,32,120,108,105,110,107,58,104,114,101,102,61,34,0,0,0,47,97,99,99,101,110,116,55,47,50,0,0,0,0,0,0,103,114,101,121,55,52,0,0,103,114,101,121,55,51,0,0,111,114,100,109,0,0,0,0,103,114,101,121,55,50,0,0,103,114,101,121,55,49,0,0,103,114,101,121,55,48,0,0,103,114,101,121,55,0,0,0,47,98,117,112,117,57,47,50,0,0,0,0,0,0,0,0,103,114,101,121,54,57,0,0,103,114,101,121,54,56,0,0,114,101,100,0,0,0,0,0,32,45,115,109,111,111,116,104,32,98,101,122,105,101,114,32,0,0,0,0,0,0,0,0,103,114,101,121,54,55,0,0,103,114,101,121,54,54,0,0,105,111,115,95,98,97,115,101,58,58,99,108,101,97,114,0,37,37,66,111,117,110,100,105,110,103,66,111,120,58,0,0,103,114,101,121,54,53,0,0,99,114,105,109,115,111,110,0,103,114,101,121,54,52,0,0,111,114,100,102,0,0,0,0,103,114,101,121,54,51,0,0,32,109,111,118,101,116,111,10,0,0,0,0,0,0,0,0,103,114,101,121,54,50,0,0,103,114,101,121,54,49,0,0,100,97,114,107,98,114,111,119,110,0,0,0,0,0,0,0,103,114,101,121,54,48,0,0,47,98,117,112,117,57,47,49,0,0,0,0,0,0,0,0,84,105,109,101,115,45,66,111,108,100,0,0,0,0,0,0,103,114,101,121,54,0,0,0,103,114,101,121,53,57,0,0,106,112,101,58,109,97,112,0,103,114,101,121,53,56,0,0,100,105,109,101,110,0,0,0,103,114,101,121,53,55,0,0,32,105,100,61,34,97,95,0,32,105,110,32,37,115,32,45,32,115,101,116,116,105,110,103,32,116,111,32,37,46,48,50,102,10,0,0,0,0,0,0,103,114,101,121,53,54,0,0,103,114,101,121,53,53,0,0,111,114,0,0,0,0,0,0,103,114,101,121,53,52,0,0,103,114,101,121,53,51,0,0,103,114,101,121,53,50,0,0,95,112,111,114,116,95,37,115,95,40,37,100,41,95,40,37,100,41,95,37,108,100,0,0,103,114,101,121,53,49,0,0,47,98,117,112,117,56,47,56,0,0,0,0,0,0,0,0,108,101,118,101,108,32,110,111,100,101,32,114,101,99,0,0,99,111,110,115,116,114,97,105,110,105,110,103,95,102,108,97,116,95,101,100,103,101,40,103,44,118,44,101,41,32,61,61,32,70,65,76,83,69,0,0,103,114,101,121,53,48,0,0,103,114,101,121,53,0,0,0,103,114,101,121,52,57,0,0,103,114,101,121,52,56,0,0,60,103,0,0,0,0,0,0,103,114,101,121,52,55,0,0,103,114,101,121,52,54,0,0,111,112,108,117,115,0,0,0,73,110,118,97,108,105,100,32,50,45,98,121,116,101,32,85,84,70,56,32,102,111,117,110,100,32,105,110,32,105,110,112,117,116,32,111,102,32,103,114,97,112,104,32,37,115,32,45,32,116,114,101,97,116,101,100,32,97,115,32,76,97,116,105,110,45,49,46,32,80,101,114,104,97,112,115,32,34,45,71,99,104,97,114,115,101,116,61,108,97,116,105,110,49,34,32,105,115,32,110,101,101,100,101,100,63,10,0,0,0,0,0,103,114,101,121,52,53,0,0,103,114,101,121,52,52,0,0,100,105,97,109,111,110,100,0,103,114,101,121,52,51,0,0,71,0,0,0,0,0,0,0,37,46,53,103,44,37,46,53,103,0,0,0,0,0,0,0,103,114,101,121,52,50,0,0,47,98,117,112,117,56,47,55,0,0,0,0,0,0,0,0,103,114,101,121,52,49,0,0,103,114,101,121,52,48,0,0,103,114,101,121,52,0,0,0,115,101,114,105,102,0,0,0,103,114,101,121,51,57,0,0,60,47,103,62,10,0,0,0,103,114,101,121,51,56,0,0,103,114,101,121,51,55,0,0,111,109,105,99,114,111,110,0,103,114,101,121,51,54,0,0,103,114,101,121,51,53,0,0,103,114,101,121,51,52,0,0,103,114,101,121,51,51,0,0,47,98,117,112,117,56,47,54,0,0,0,0,0,0,0,0,103,114,101,121,51,50,0,0,103,114,101,121,51,49,0,0,103,114,101,121,51,48,0,0,103,114,101,121,51,0,0,0,37,37,37,37,67,114,101,97,116,111,114,58,32,37,115,32,118,101,114,115,105,111,110,32,37,115,32,40,37,115,41,10,0,0,0,0,0,0,0,0,103,114,101,121,50,57,0,0,103,114,101,121,50,56,0,0,111,109,101,103,97,0,0,0,103,114,101,121,50,55,0,0,103,114,101,121,50,54,0,0,103,114,101,121,50,53,0,0,103,114,101,121,50,52,0,0,47,98,117,112,117,56,47,53,0,0,0,0,0,0,0,0,103,114,101,121,50,51,0,0,103,114,101,121,50,50,0,0,103,114,101,121,50,49,0,0,103,114,101,121,50,48,0,0,60,47,116,101,120,116,62,10,0,0,0,0,0,0,0,0,32,69,80,83,70,45,51,46,48,10,0,0,0,0,0,0,103,114,101,121,50,0,0,0,103,114,101,121,49,57,0,0,111,108,105,110,101,0,0,0,103,114,101,121,49,56,0,0,103,114,101,121,49,55,0,0,103,114,101,121,49,54,0,0,103,114,101,121,49,53,0,0,47,98,117,112,117,56,47,52,0,0,0,0,0,0,0,0,97,103,110,111,100,101,0,0,103,114,101,121,49,52,0,0,103,114,101,121,49,51,0,0,103,114,101,121,49,50,0,0,103,114,101,121,49,49,0,0,62,0,0,0,0,0,0,0,37,33,80,83,45,65,100,111,98,101,45,51,46,48,0,0,103,114,101,121,49,48,48,0,103,114,101,121,49,48,0,0,111,103,114,97,118,101,0,0,103,114,101,121,49,0,0,0,103,114,101,121,48,0,0,0,47,98,117,112,117,56,47,51,0,0,0,0,0,0,0,0,103,114,101,101,110,52,0,0,103,114,101,101,110,51,0,0,103,114,101,101,110,50,0,0,103,114,101,101,110,49,0,0,32,102,105,108,108,61,34,35,37,48,50,120,37,48,50,120,37,48,50,120,34,0,0,0,37,37,69,79,70,10,0,0,103,114,97,121,57,57,0,0,111,101,108,105,103,0,0,0,103,114,97,121,57,56,0,0,103,114,97,121,57,55,0,0,103,114,97,121,57,54,0,0,47,98,117,112,117,56,47,50,0,0,0,0,0,0,0,0,103,114,97,121,57,52,0,0,103,114,97,121,57,51,0,0,103,114,97,121,57,50,0,0,103,114,97,121,57,49,0,0,32,102,105,108,108,61,34,37,115,34,0,0,0,0,0,0,101,110,100,10,114,101,115,116,111,114,101,10,0,0,0,0,103,114,97,121,57,0,0,0,111,99,105,114,99,0,0,0,103,114,97,121,56,57,0,0,103,114,97,121,56,56,0,0,103,114,97,121,56,55,0,0,103,114,97,121,56,54,0,0,47,98,117,112,117,56,47,49,0,0,0,0,0,0,0,0,103,114,97,121,56,52,0,0,103,114,97,121,56,51,0,0,103,114,97,121,56,50,0,0,32,102,111,110,116,45,115,105,122,101,61,34,37,46,50,102,34,0,0,0,0,0,0,0,37,37,37,37,80,97,103,101,115,58,32,37,100,10,0,0,47,97,99,99,101,110,116,55,47,49,0,0,0,0,0,0,103,114,97,121,56,49,0,0,111,97,99,117,116,101,0,0,103,114,97,121,56,0,0,0,103,114,97,121,55,57,0,0,103,114,97,121,55,56,0,0,103,114,97,121,55,55,0,0,47,98,117,112,117,55,47,55,0,0,0,0,0,0,0,0,103,114,97,121,55,54,0,0,112,117,114,112,108,101,0,0,32,45,119,105,100,116,104,32,0,0,0,0,0,0,0,0,103,114,97,121,55,52,0,0,103,114,97,121,55,51,0,0,32,98,97,115,101,108,105,110,101,45,115,104,105,102,116,61,34,115,117,98,34,0,0,0,37,37,84,114,97,105,108,101,114,10,0,0,0,0,0,0,40,91,97,45,122,93,91,97,45,122,65,45,90,93,42,41,61,34,40,91,94,34,93,42,41,34,0,0,0,0,0,0,103,114,97,121,55,50,0,0,99,111,114,110,115,105,108,107,0,0,0,0,0,0,0,0,103,114,97,121,55,49,0,0,110,117,0,0,0,0,0,0,110,101,119,112,97,116,104,32,0,0,0,0,0,0,0,0,103,114,97,121,55,0,0,0,103,114,97,121,54,57,0,0,103,114,97,121,54,56,0,0,47,98,117,112,117,55,47,54,0,0,0,0,0,0,0,0,65,118,97,110,116,71,97,114,100,101,45,68,101,109,105,79,98,108,105,113,117,101,0,0,32,104,114,101,102,61,34,0,37,100,32,37,100,32,35,37,48,50,120,37,48,50,120,37,48,50,120,10,0,0,0,0,103,114,97,121,54,55,0,0,103,114,97,121,54,54,0,0,106,112,101,103,58,109,97,112,0,0,0,0,0,0,0,0,115,112,114,105,110,103,0,0,103,114,97,121,54,52,0,0,32,98,97,115,101,108,105,110,101,45,115,104,105,102,116,61,34,115,117,112,101,114,34,0,37,37,69,110,100,83,101,116,117,112,0,0,0,0,0,0,103,114,97,121,54,51,0,0,103,114,97,121,54,50,0,0,110,116,105,108,100,101,0,0,103,114,97,112,104,32,37,115,32,105,115,32,100,105,115,99,111,110,110,101,99,116,101,100,46,32,72,101,110,99,101,44,32,116,104,101,32,99,105,114,99,117,105,116,32,109,111,100,101,108,10,0,0,0,0,0,103,114,97,121,54,49,0,0,103,114,97,121,54,0,0,0,111,118,101,114,108,97,112,95,115,99,97,108,105,110,103,0,95,112,111,114,116,95,37,115,95,37,115,95,37,115,95,37,108,100,0,0,0,0,0,0,103,114,97,121,53,57,0,0,47,98,117,112,117,55,47,53,0,0,0,0,0,0,0,0,108,101,118,101,108,32,101,100,103,101,32,114,101,99,0,0,103,114,97,121,53,56,0,0,103,114,97,121,53,55,0,0,103,114,97,121,53,54,0,0,32,116,101,120,116,45,100,101,99,111,114,97,116,105,111,110,61,34,117,110,100,101,114,108,105,110,101,34,0,0,0,0,125,32,105,102,0,0,0,0,103,114,97,121,53,52,0,0,103,114,97,121,53,51,0,0,110,115,117,98,0,0,0,0,103,114,97,121,53,50,0,0,103,114,97,121,53,49,0,0,112,108,97,105,110,116,101,120,116,0,0,0,0,0,0,0,108,104,101,105,103,104,116,0,70,0,0,0,0,0,0,0,92,76,0,0,0,0,0,0,103,114,97,121,53,0,0,0,47,98,117,112,117,55,47,52,0,0,0,0,0,0,0,0,123,10,0,0,0,0,0,0,103,114,97,121,52,57,0,0,103,114,97,121,52,56,0,0,103,114,97,121,52,55,0,0,85,82,87,32,66,111,111,107,109,97,110,32,76,0,0,0,103,114,97,121,52,54,0,0,32,102,111,110,116,45,115,116,121,108,101,61,34,105,116,97,108,105,99,34,0,0,0,0,32,32,32,32,117,115,101,114,100,105,99,116,32,40,62,62,41,32,99,118,110,32,40,91,41,32,99,118,110,32,108,111,97,100,32,112,117,116,0,0,103,114,97,121,52,52,0,0,110,111,116,105,110,0,0,0,103,114,97,121,52,51,0,0,103,114,97,121,52,50,0,0,103,114,97,121,52,49,0,0,47,98,117,112,117,55,47,51,0,0,0,0,0,0,0,0,103,114,97,121,52,0,0,0,103,114,97,121,51,57,0,0,103,114,97,121,51,56,0,0,103,114,97,121,51,55,0,0,32,102,111,110,116,45,119,101,105,103,104,116,61,34,98,111,108,100,34,0,0,0,0,0,32,32,32,32,117,115,101,114,100,105,99,116,32,40,60,60,41,32,99,118,110,32,40,91,41,32,99,118,110,32,108,111,97,100,32,112,117,116,0,0,103,114,97,121,51,54,0,0,110,111,116,0,0,0,0,0,103,114,97,121,51,52,0,0,103,114,97,121,51,51,0,0,103,114,97,121,51,50,0,0,103,114,97,121,51,49,0,0,47,98,117,112,117,55,47,50,0,0,0,0,0,0,0,0,103,114,97,121,51,0,0,0,103,114,97,121,50,57,0,0,103,114,97,121,50,56,0,0,32,102,111,110,116,45,102,97,109,105,108,121,61,34,37,115,34,0,0,0,0,0,0,0,50,32,108,116,32,123,0,0,103,114,97,121,50,55,0,0,103,114,97,121,50,54,0,0,110,105,0,0,0,0,0,0,103,114,97,121,50,52,0,0,103,114,97,121,50,51,0,0,103,114,97,121,50,50,0,0,47,98,117,112,117,55,47,49,0,0,0,0,0,0,0,0,97,103,110,101,100,103,101,115,0,0,0,0,0,0,0,0,103,114,97,121,50,49,0,0,47,98,117,112,117,54,47,54,0,0,0,0,0,0,0,0,103,114,97,121,50,0,0,0,103,114,97,121,49,57,0,0,32,102,111,110,116,45,115,116,121,108,101,61,34,37,115,34,0,0,0,0,0,0,0,0,47,108,97,110,103,117,97,103,101,108,101,118,101,108,32,119,104,101,114,101,32,123,112,111,112,32,108,97,110,103,117,97,103,101,108,101,118,101,108,125,123,49,125,32,105,102,101,108,115,101,0,0,0,0,0,0,103,114,97,121,49,56,0,0,103,114,97,121,49,55,0,0,110,101,0,0,0,0,0,0,103,114,97,121,49,54,0,0,103,114,97,121,49,52,0,0,103,114,97,121,49,51,0,0,109,111,118,101,32,116,111,32,102,114,111,110,116,32,108,111,99,107,32,105,110,99,111,110,115,105,115,116,101,110,99,121,0,0,0,0,0,0,0,0,103,114,97,121,49,50,0,0,103,114,97,121,49,49,0,0,103,114,97,121,49,48,48,0,32,102,111,110,116,45,115,116,114,101,116,99,104,61,34,37,115,34,0,0,0,0,0,0,37,32,109,97,107,101,32,39,60,60,39,32,97,110,100,32,39,62,62,39,32,115,97,102,101,32,111,110,32,80,83,32,76,101,118,101,108,32,49,32,100,101,118,105,99,101,115,0,103,114,97,121,49,0,0,0,103,114,97,121,48,0,0,0,110,100,97,115,104,0,0,0,103,111,108,100,101,110,114,111,100,52,0,0,0,0,0,0,103,111,108,100,101,110,114,111,100,51,0,0,0,0,0,0,103,111,108,100,101,110,114,111,100,50,0,0,0,0,0,0,47,98,117,112,117,54,47,53,0,0,0,0,0,0,0,0,103,111,108,100,101,110,114,111,100,49,0,0,0,0,0,0,103,111,108,100,52,0,0,0,103,111,108,100,51,0,0,0,32,102,111,110,116,45,119,101,105,103,104,116,61,34,37,115,34,0,0,0,0,0,0,0,47,112,100,102,109,97,114,107,32,119,104,101,114,101,32,123,112,111,112,125,32,123,117,115,101,114,100,105,99,116,32,47,112,100,102,109,97,114,107,32,47,99,108,101,97,114,116,111,109,97,114,107,32,108,111,97,100,32,112,117,116,125,32,105,102,101,108,115,101,0,0,0,103,111,108,100,50,0,0,0,103,111,108,100,49,0,0,0,110,98,115,112,0,0,0,0,47,98,117,112,117,54,47,52,0,0,0,0,0,0,0,0,102,105,114,101,98,114,105,99,107,52,0,0,0,0,0,0,102,105,114,101,98,114,105,99,107,51,0,0,0,0,0,0,102,105,114,101,98,114,105,99,107,50,0,0,0,0,0,0,44,37,115,0,0,0,0,0,37,32,109,97,107,101,32,115,117,114,101,32,112,100,102,109,97,114,107,32,105,115,32,104,97,114,109,108,101,115,115,32,102,111,114,32,80,83,45,105,110,116,101,114,112,114,101,116,101,114,115,32,111,116,104,101,114,32,116,104,97,110,32,68,105,115,116,105,108,108,101,114,0,0,0,0,0,0,0,0,47,97,99,99,101,110,116,54,47,54,0,0,0,0,0,0,102,105,114,101,98,114,105,99,107,49,0,0,0,0,0,0,110,97,98,108,97,0,0,0,100,111,100,103,101,114,98,108,117,101,52,0,0,0,0,0,100,111,100,103,101,114,98,108,117,101,51,0,0,0,0,0,100,111,100,103,101,114,98,108,117,101,50,0,0,0,0,0,100,111,100,103,101,114,98,108,117,101,49,0,0,0,0,0,47,98,117,112,117,54,47,51,0,0,0,0,0,0,0,0,37,108,100,0,0,0,0,0,111,108,105,118,101,0,0,0,36,99,0,0,0,0,0,0,100,101,101,112,115,107,121,98,108,117,101,52,0,0,0,0,32,102,111,110,116,45,102,97,109,105,108,121,61,34,37,115,0,0,0,0,0,0,0,0,37,32,47,97,114,114,111,119,119,105,100,116,104,32,53,32,100,101,102,0,0,0,0,0,100,101,101,112,115,107,121,98,108,117,101,51,0,0,0,0,109,109,0,0,0,0,0,0,99,111,114,110,102,108,111,119,101,114,98,108,117,101,0,0,100,101,101,112,115,107,121,98,108,117,101,50,0,0,0,0,109,117,0,0,0,0,0,0,100,101,101,112,115,107,121,98,108,117,101,49,0,0,0,0,100,101,101,112,112,105,110,107,52,0,0,0,0,0,0,0,100,101,101,112,112,105,110,107,51,0,0,0,0,0,0,0,47,98,117,112,117,54,47,50,0,0,0,0,0,0,0,0,65,118,97,110,116,71,97,114,100,101,45,66,111,111,107,0,37,100,32,37,100,32,37,100,32,37,100,32,37,100,32,37,100,32,37,100,32,37,100,32,37,100,32,37,46,51,102,32,37,100,32,37,46,52,102,32,37,100,32,37,100,32,37,100,32,37,100,32,37,100,32,37,100,32,37,100,32,37,100,10,0,0,0,0,0,0,0,0,100,101,101,112,112,105,110,107,50,0,0,0,0,0,0,0,67,32,0,0,0,0,0,0,100,101,101,112,112,105,110,107,49,0,0,0,0,0,0,0,103,105,102,58,109,97,112,0,47,98,117,112,117,54,47,49,0,0,0,0,0,0,0,0])
.concat([112,111,119,101,114,95,100,105,115,116,0,0,0,0,0,0,32,120,61,34,37,103,34,32,121,61,34,37,103,34,0,0,37,32,47,97,114,114,111,119,108,101,110,103,116,104,32,49,48,32,100,101,102,0,0,0,78,68,95,104,101,97,112,105,110,100,101,120,40,118,41,32,60,32,48,0,0,0,0,0,109,105,110,117,115,0,0,0,109,97,120,105,116,101,114,0,100,97,114,107,115,108,97,116,101,103,114,97,121,52,0,0,100,97,114,107,115,108,97,116,101,103,114,97,121,51,0,0,100,97,114,107,115,108,97,116,101,103,114,97,121,50,0,0,110,111,100,101,32,34,37,115,34,32,105,115,32,99,111,110,116,97,105,110,101,100,32,105,110,32,116,119,111,32,110,111,110,45,99,111,109,112,97,114,97,98,108,101,32,99,108,117,115,116,101,114,115,32,34,37,115,34,32,97,110,100,32,34,37,115,34,10,0,0,0,0,100,97,114,107,115,108,97,116,101,103,114,97,121,49,0,0,97,103,114,101,99,111,114,100,95,99,97,108,108,98,97,99,107,32,111,102,32,97,32,98,97,100,32,111,98,106,101,99,116,0,0,0,0,0,0,0,115,97,109,101,0,0,0,0,100,97,114,107,115,101,97,103,114,101,101,110,52,0,0,0,100,97,114,107,115,101,97,103,114,101,101,110,51,0,0,0,32,116,101,120,116,45,97,110,99,104,111,114,61,34,109,105,100,100,108,101,34,0,0,0,49,32,115,101,116,109,105,116,101,114,108,105,109,105,116,0,100,97,114,107,115,101,97,103,114,101,101,110,50,0,0,0,100,97,114,107,115,101,97,103,114,101,101,110,49,0,0,0,109,105,100,100,111,116,0,0,101,110,100,32,112,111,114,116,58,32,40,37,46,53,103,44,32,37,46,53,103,41,44,32,116,97,110,103,101,110,116,32,97,110,103,108,101,58,32,37,46,53,103,44,32,37,115,10,0,0,0,0,0,0,0,0,100,97,114,107,111,114,99,104,105,100,52,0,0,0,0,0,108,119,105,100,116,104,0,0,92,84,0,0,0,0,0,0,100,97,114,107,111,114,99,104,105,100,51,0,0,0,0,0,47,98,117,112,117,53,47,53,0,0,0,0,0,0,0,0,100,97,114,107,111,114,99,104,105,100,50,0,0,0,0,0,100,97,114,107,111,114,99,104,105,100,49,0,0,0,0,0,100,97,114,107,111,114,97,110,103,101,52,0,0,0,0,0,32,116,101,120,116,45,97,110,99,104,111,114,61,34,101,110,100,34,0,0,0,0,0,0,49,52,32,100,101,102,97,117,108,116,45,102,111,110,116,45,102,97,109,105,108,121,32,115,101,116,95,102,111,110,116,0,100,97,114,107,111,114,97,110,103,101,51,0,0,0,0,0,37,100,32,37,49,91,34,93,37,110,0,0,0,0,0,0,100,97,114,107,111,114,97,110,103,101,50,0,0,0,0,0,109,105,99,114,111,0,0,0,100,97,114,107,111,114,97,110,103,101,49,0,0,0,0,0,100,97,114,107,111,108,105,118,101,103,114,101,101,110,52,0,100,97,114,107,111,108,105,118,101,103,114,101,101,110,51,0,47,98,117,112,117,53,47,52,0,0,0,0,0,0,0,0,100,97,114,107,111,108,105,118,101,103,114,101,101,110,50,0,100,97,114,107,111,108,105,118,101,103,114,101,101,110,49,0,32,116,101,120,116,45,97,110,99,104,111,114,61,34,115,116,97,114,116,34,0,0,0,0,37,37,66,101,103,105,110,83,101,116,117,112,0,0,0,0,100,97,114,107,103,111,108,100,101,110,114,111,100,52,0,0,109,100,97,115,104,0,0,0,100,97,114,107,103,111,108,100,101,110,114,111,100,51,0,0,100,97,114,107,103,111,108,100,101,110,114,111,100,50,0,0,100,97,114,107,103,111,108,100,101,110,114,111,100,49,0,0,47,98,117,112,117,53,47,51,0,0,0,0,0,0,0,0,99,121,97,110,52,0,0,0,99,121,97,110,51,0,0,0,99,121,97,110,50,0,0,0,99,121,97,110,49,0,0,0,60,116,101,120,116,0,0,0,37,37,69,110,100,80,114,111,108,111,103,0,0,0,0,0,109,97,99,114,0,0,0,0,99,111,114,110,115,105,108,107,52,0,0,0,0,0,0,0,99,111,114,110,115,105,108,107,51,0,0,0,0,0,0,0,99,111,114,110,115,105,108,107,50,0,0,0,0,0,0,0,99,111,114,110,115,105,108,107,49,0,0,0,0,0,0,0,47,98,117,112,117,53,47,50,0,0,0,0,0,0,0,0,97,103,110,110,111,100,101,115,0,0,0,0,0,0,0,0,99,111,114,97,108,52,0,0,99,111,114,97,108,51,0,0,47,62,10,0,0,0,0,0,37,37,69,110,100,82,101,115,111,117,114,99,101,0,0,0,99,111,114,97,108,50,0,0,99,111,114,97,108,49,0,0,108,116,0,0,0,0,0,0,99,104,111,99,111,108,97,116,101,52,0,0,0,0,0,0,99,104,111,99,111,108,97,116,101,51,0,0,0,0,0,0,99,104,111,99,111,108,97,116,101,50,0,0,0,0,0,0,47,98,117,112,117,53,47,49,0,0,0,0,0,0,0,0,99,104,111,99,111,108,97,116,101,49,0,0,0,0,0,0,99,104,97,114,116,114,101,117,115,101,52,0,0,0,0,0,99,104,97,114,116,114,101,117,115,101,51,0,0,0,0,0,32,114,120,61,34,37,103,34,32,114,121,61,34,37,103,34,0,0,0,0,0,0,0,0,47,99,117,114,108,97,121,101,114,32,48,32,100,101,102,0,99,104,97,114,116,114,101,117,115,101,50,0,0,0,0,0,99,104,97,114,116,114,101,117,115,101,49,0,0,0,0,0,108,115,113,117,111,0,0,0,99,97,100,101,116,98,108,117,101,52,0,0,0,0,0,0,99,97,100,101,116,98,108,117,101,51,0,0,0,0,0,0,99,97,100,101,116,98,108,117,101,50,0,0,0,0,0,0,47,98,117,112,117,52,47,52,0,0,0,0,0,0,0,0,99,97,100,101,116,98,108,117,101,49,0,0,0,0,0,0,98,117,114,108,121,119,111,111,100,52,0,0,0,0,0,0,47,98,117,112,117,52,47,51,0,0,0,0,0,0,0,0,98,117,114,108,121,119,111,111,100,51,0,0,0,0,0,0,32,99,120,61,34,37,103,34,32,99,121,61,34,37,103,34,0,0,0,0,0,0,0,0,9,123,105,110,118,105,115,125,32,105,102,0,0,0,0,0,98,117,114,108,121,119,111,111,100,50,0,0,0,0,0,0,98,117,114,108,121,119,111,111,100,49,0,0,0,0,0,0,108,115,97,113,117,111,0,0,98,114,111,119,110,52,0,0,98,114,111,119,110,51,0,0,98,114,111,119,110,50,0,0,98,114,111,119,110,49,0,0,98,108,117,101,52,0,0,0,60,101,108,108,105,112,115,101,0,0,0,0,0,0,0,0,9,111,114,0,0,0,0,0,47,97,99,99,101,110,116,54,47,53,0,0,0,0,0,0,98,108,117,101,51,0,0,0,98,108,117,101,50,0,0,0,108,114,109,0,0,0,0,0,98,108,117,101,49,0,0,0,47,98,117,112,117,52,47,50,0,0,0,0,0,0,0,0,98,105,115,113,117,101,52,0,98,105,115,113,117,101,51,0,110,97,118,121,0,0,0,0,98,105,115,113,117,101,50,0,98,105,115,113,117,101,49,0,37,103,44,37,103,0,0,0,9,99,117,114,108,97,121,101,114,32,109,121,117,112,112,101,114,32,103,116,0,0,0,0,99,111,114,97,108,0,0,0,99,109,0,0,0,0,0,0,108,111,122,0,0,0,0,0,97,122,117,114,101,52,0,0,37,32,0,0,0,0,0,0,97,122,117,114,101,51,0,0,97,122,117,114,101,50,0,0,97,122,117,114,101,49,0,0,47,98,117,112,117,52,47,49,0,0,0,0,0,0,0,0,65,118,97,110,116,71,97,114,100,101,45,66,111,111,107,79,98,108,105,113,117,101,0,0,32,105,100,61,34,0,0,0,97,113,117,97,109,97,114,105,110,101,52,0,0,0,0,0,112,110,103,58,109,97,112,0,97,113,117,97,109,97,114,105,110,101,51,0,0,0,0,0,103,114,97,112,104,95,100,105,115,116,0,0,0,0,0,0,97,113,117,97,109,97,114,105,110,101,50,0,0,0,0,0,60,112,111,108,121,103,111,110,0,0,0,0,0,0,0,0,9,99,117,114,108,97,121,101,114,32,109,121,108,111,119,101,114,32,108,116,0,0,0,0,46,46,47,46,46,47,46,46,47,108,105,98,47,110,101,97,116,111,103,101,110,47,115,116,117,102,102,46,99,0,0,0,97,113,117,97,109,97,114,105,110,101,49,0,0,0,0,0,108,111,119,97,115,116,0,0,97,110,116,105,113,117,101,119,104,105,116,101,52,0,0,0,97,110,116,105,113,117,101,119,104,105,116,101,51,0,0,0,97,110,116,105,113,117,101,119,104,105,116,101,50,0,0,0,79,118,101,114,108,97,112,32,118,97,108,117,101,32,34,37,115,34,32,117,110,115,117,112,112,111,114,116,101,100,32,45,32,105,103,110,111,114,101,100,10,0,0,0,0,0,0,0,97,110,116,105,113,117,101,119,104,105,116,101,49,0,0,0,47,98,117,112,117,51,47,51,0,0,0,0,0,0,0,0,115,105,110,107,0,0,0,0,40,108,32,61,32,69,68,95,108,97,98,101,108,40,102,101,41,41,0,0,0,0,0,0,47,121,108,111,114,114,100,57,47,57,0,0,0,0,0,0,47,121,108,111,114,114,100,57,47,56,0,0,0,0,0,0,59,34,47,62,10,60,47,108,105,110,101,97,114,71,114,97,100,105,101,110,116,62,10,60,47,100,101,102,115,62,10,0,9,47,109,121,108,111,119,101,114,32,101,120,99,104,32,100,101,102,0,0,0,0,0,0,87,104,105,116,101,0,0,0,47,121,108,111,114,114,100,57,47,55,0,0,0,0,0,0,47,121,108,111,114,114,100,57,47,54,0,0,0,0,0,0,108,102,108,111,111,114,0,0,95,95,99,108,117,115,116,101,114,110,111,100,101,115,0,0,47,121,108,111,114,114,100,57,47,53,0,0,0,0,0,0,47,121,108,111,114,114,100,57,47,52,0,0,0,0,0,0,116,114,105,97,110,103,108,101,0,0,0,0,0,0,0,0,110,111,116,32,99,111,110,115,116,114,97,105,110,101,100,0,66,111,117,110,100,105,110,103,66,111,120,32,110,111,116,32,102,111,117,110,100,32,105,110,32,101,112,115,102,32,102,105,108,101,32,37,115,10,0,0,47,121,108,111,114,114,100,57,47,51,0,0,0,0,0,0,92,72,0,0,0,0,0,0,47,121,108,111,114,114,100,57,47,50,0,0,0,0,0,0,47,98,117,112,117,51,47,50,0,0,0,0,0,0,0,0,47,121,108,111,114,114,100,57,47,49,0,0,0,0,0,0,47,121,108,111,114,114,100,56,47,56,0,0,0,0,0,0,47,121,108,111,114,114,100,56,47,55,0,0,0,0,0,0,47,121,108,111,114,114,100,56,47,54,0,0,0,0,0,0,120,49,61,34,37,103,34,32,121,49,61,34,37,103,34,32,120,50,61,34,37,103,34,32,121,50,61,34,37,103,34,32,62,10,0,0,0,0,0,0,9,47,109,121,117,112,112,101,114,32,101,120,99,104,32,100,101,102,0,0,0,0,0,0,66,108,97,99,107,0,0,0,47,121,108,111,114,114,100,56,47,53,0,0,0,0,0,0,108,105,110,101,0,0,0,0,47,121,108,111,114,114,100,56,47,52,0,0,0,0,0,0,108,101,0,0,0,0,0,0,47,121,108,111,114,114,100,56,47,51,0,0,0,0,0,0,47,121,108,111,114,114,100,56,47,50,0,0,0,0,0,0,47,121,108,111,114,114,100,56,47,49,0,0,0,0,0,0,47,121,108,111,114,114,100,55,47,55,0,0,0,0,0,0,47,98,117,112,117,51,47,49,0,0,0,0,0,0,0,0,47,121,108,111,114,114,100,55,47,54,0,0,0,0,0,0,47,121,108,111,114,114,100,55,47,53,0,0,0,0,0,0,47,121,108,111,114,114,100,55,47,52,0,0,0,0,0,0,47,121,108,111,114,114,100,55,47,51,0,0,0,0,0,0,60,100,101,102,115,62,10,60,108,105,110,101,97,114,71,114,97,100,105,101,110,116,32,105,100,61,34,108,95,37,100,34,32,103,114,97,100,105,101,110,116,85,110,105,116,115,61,34,117,115,101,114,83,112,97,99,101,79,110,85,115,101,34,32,0,0,0,0,0,0,0,0,47,111,110,108,97,121,101,114,115,32,123,0,0,0,0,0,35,100,101,99,108,97,114,101,32,37,115,32,61,32,37,115,59,10,0,0,0,0,0,0,47,121,108,111,114,114,100,55,47,50,0,0,0,0,0,0,47,121,108,111,114,114,100,55,47,49,0,0,0,0,0,0,108,100,113,117,111,0,0,0,47,121,108,111,114,114,100,54,47,54,0,0,0,0,0,0,47,121,108,111,114,114,100,54,47,53,0,0,0,0,0,0,47,121,108,111,114,114,100,54,47,52,0,0,0,0,0,0,47,121,108,111,114,114,100,54,47,51,0,0,0,0,0,0,47,98,117,103,110,57,47,57,0,0,0,0,0,0,0,0,47,121,108,111,114,114,100,54,47,50,0,0,0,0,0,0,47,121,108,111,114,114,100,54,47,49,0,0,0,0,0,0,47,121,108,111,114,114,100,53,47,53,0,0,0,0,0,0,47,121,108,111,114,114,100,53,47,52,0,0,0,0,0,0,59,34,47,62,10,60,47,114,97,100,105,97,108,71,114,97,100,105,101,110,116,62,10,60,47,100,101,102,115,62,10,0,47,111,110,108,97,121,101,114,32,123,32,99,117,114,108,97,121,101,114,32,110,101,32,123,105,110,118,105,115,125,32,105,102,32,125,32,100,101,102,0,35,105,110,99,108,117,100,101,32,34,99,111,108,111,114,115,46,105,110,99,34,10,35,105,110,99,108,117,100,101,32,34,116,101,120,116,117,114,101,115,46,105,110,99,34,10,35,105,110,99,108,117,100,101,32,34,115,104,97,112,101,115,46,105,110,99,34,10,0,0,0,0,47,121,108,111,114,114,100,53,47,51,0,0,0,0,0,0,47,121,108,111,114,114,100,53,47,50,0,0,0,0,0,0,108,99,101,105,108,0,0,0,47,121,108,111,114,114,100,53,47,49,0,0,0,0,0,0,47,121,108,111,114,114,100,52,47,52,0,0,0,0,0,0,47,121,108,111,114,114,100,52,47,51,0,0,0,0,0,0,47,121,108,111,114,114,100,52,47,50,0,0,0,0,0,0,47,98,117,103,110,57,47,56,0,0,0,0,0,0,0,0,97,103,99,111,110,99,97,116,0,0,0,0,0,0,0,0,47,121,108,111,114,114,100,52,47,49,0,0,0,0,0,0,47,121,108,111,114,114,100,51,47,51,0,0,0,0,0,0,47,121,108,111,114,114,100,51,47,50,0,0,0,0,0,0,47,121,108,111,114,114,100,51,47,49,0,0,0,0,0,0,60,115,116,111,112,32,111,102,102,115,101,116,61,34,49,34,32,115,116,121,108,101,61,34,115,116,111,112,45,99,111,108,111,114,58,0,0,0,0,0,9,47,103,114,97,112,104,99,111,108,111,114,32,123,110,111,112,99,111,108,111,114,125,32,100,101,102,0,0,0,0,0,35,100,101,102,97,117,108,116,32,123,32,102,105,110,105,115,104,32,123,32,97,109,98,105,101,110,116,32,48,46,49,32,100,105,102,102,117,115,101,32,48,46,57,32,125,32,125,10,0,0,0,0,0,0,0,0,47,121,108,111,114,98,114,57,47,57,0,0,0,0,0,0,47,121,108,111,114,98,114,57,47,56,0,0,0,0,0,0,108,97,114,114,0,0,0,0,47,121,108,111,114,98,114,57,47,55,0,0,0,0,0,0,47,121,108,111,114,98,114,57,47,54,0,0,0,0,0,0,47,121,108,111,114,98,114,57,47,53,0,0,0,0,0,0,47,121,108,111,114,98,114,57,47,52,0,0,0,0,0,0,47,98,117,103,110,57,47,55,0,0,0,0,0,0,0,0,47,121,108,111,114,98,114,57,47,51,0,0,0,0,0,0,47,121,108,111,114,98,114,57,47,50,0,0,0,0,0,0,47,121,108,111,114,98,114,57,47,49,0,0,0,0,0,0,47,98,117,103,110,57,47,54,0,0,0,0,0,0,0,0,47,121,108,111,114,98,114,56,47,56,0,0,0,0,0,0,59,34,47,62,10,0,0,0,9,47,101,100,103,101,99,111,108,111,114,32,123,110,111,112,99,111,108,111,114,125,32,100,101,102,0,0,0,0,0,0,103,108,111,98,97,108,95,115,101,116,116,105,110,103,115,32,123,32,97,115,115,117,109,101,100,95,103,97,109,109,97,32,49,46,48,32,125,10,0,0,47,121,108,111,114,98,114,56,47,55,0,0,0,0,0,0,47,121,108,111,114,98,114,56,47,54,0,0,0,0,0,0,108,97,113,117,111,0,0,0,47,121,108,111,114,98,114,56,47,53,0,0,0,0,0,0,47,121,108,111,114,98,114,56,47,52,0,0,0,0,0,0,47,121,108,111,114,98,114,56,47,51,0,0,0,0,0,0,47,121,108,111,114,98,114,56,47,50,0,0,0,0,0,0,97,103,100,101,108,101,116,101,32,111,110,32,119,114,111,110,103,32,103,114,97,112,104,0,47,121,108,111,114,98,114,56,47,49,0,0,0,0,0,0,47,121,108,111,114,98,114,55,47,55,0,0,0,0,0,0,47,121,108,111,114,98,114,55,47,54,0,0,0,0,0,0,47,121,108,111,114,98,114,55,47,53,0,0,0,0,0,0,49,46,0,0,0,0,0,0,9,47,110,111,100,101,99,111,108,111,114,32,123,110,111,112,99,111,108,111,114,125,32,100,101,102,0,0,0,0,0,0,35,118,101,114,115,105,111,110,32,51,46,54,59,10,0,0,47,121,108,111,114,98,114,55,47,52,0,0,0,0,0,0,47,121,108,111,114,98,114,55,47,51,0,0,0,0,0,0,108,97,110,103,0,0,0,0,47,121,108,111,114,98,114,55,47,50,0,0,0,0,0,0,47,121,108,111,114,98,114,55,47,49,0,0,0,0,0,0,47,121,108,111,114,98,114,54,47,54,0,0,0,0,0,0,47,121,108,111,114,98,114,54,47,53,0,0,0,0,0,0,47,98,117,103,110,57,47,53,0,0,0,0,0,0,0,0,47,121,108,111,114,98,114,54,47,52,0,0,0,0,0,0,47,121,108,111,114,98,114,54,47,51,0,0,0,0,0,0,47,121,108,111,114,98,114,54,47,50,0,0,0,0,0,0,47,121,108,111,114,98,114,54,47,49,0,0,0,0,0,0,9,97,108,111,97,100,32,112,111,112,32,115,101,116,104,115,98,99,111,108,111,114,0,0,108,105,103,104,116,95,115,111,117,114,99,101,32,123,32,60,49,53,48,48,44,51,48,48,48,44,45,50,53,48,48,62,32,99,111,108,111,114,32,87,104,105,116,101,32,125,10,0,47,121,108,111,114,98,114,53,47,53,0,0,0,0,0,0,47,97,99,99,101,110,116,54,47,52,0,0,0,0,0,0,47,121,108,111,114,98,114,53,47,52,0,0,0,0,0,0,108,97,109,98,100,97,0,0,47,121,108,111,114,98,114,53,47,51,0,0,0,0,0,0,105,110,118,101,109,112,116,121,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,47,121,108,111,114,98,114,53,47,50,0,0,0,0,0,0,47,121,108,111,114,98,114,53,47,49,0,0,0,0,0,0,47,121,108,111,114,98,114,52,47,52,0,0,0,0,0,0,47,98,117,103,110,57,47,52,0,0,0,0,0,0,0,0,37,100,0,0,0,0,0,0,47,121,108,111,114,98,114,52,47,51,0,0,0,0,0,0,47,121,108,111,114,98,114,52,47,50,0,0,0,0,0,0,109,97,114,111,111,110,0,0,34,34,0,0,0,0,0,0,47,121,108,111,114,98,114,52,47,49,0,0,0,0,0,0,47,121,108,111,114,98,114,51,47,51,0,0,0,0,0,0,59,115,116,111,112,45,111,112,97,99,105,116,121,58,0,0,9,108,97,121,101,114,99,111,108,111,114,115,101,113,32,99,117,114,108,97,121,101,114,32,49,32,115,117,98,32,108,97,121,101,114,108,101,110,32,109,111,100,32,103,101,116,0,0,47,47,115,107,121,10,112,108,97,110,101,32,123,32,60,48,44,32,49,44,32,48,62,44,32,49,32,104,111,108,108,111,119,10,32,32,32,32,116,101,120,116,117,114,101,32,123,10,32,32,32,32,32,32,32,32,112,105,103,109,101,110,116,32,123,32,98,111,122,111,32,116,117,114,98,117,108,101,110,99,101,32,48,46,57,53,10,32,32,32,32,32,32,32,32,32,32,32,32,99,111,108,111,114,95,109,97,112,32,123,10,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,91,48,46,48,48,32,114,103,98,32,60,48,46,48,53,44,32,48,46,50,48,44,32,48,46,53,48,62,93,10,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,91,48,46,53,48,32,114,103,98,32,60,48,46,48,53,44,32,48,46,50,48,44,32,48,46,53,48,62,93,10,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,91,48,46,55,53,32,114,103,98,32,60,49,46,48,48,44,32,49,46,48,48,44,32,49,46,48,48,62,93,10,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,91,48,46,55,53,32,114,103,98,32,60,48,46,50,53,44,32,48,46,50,53,44,32,48,46,50,53,62,93,10,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,91,49,46,48,48,32,114,103,98,32,60,48,46,53,48,44,32,48,46,53,48,44,32,48,46,53,48,62,93,10,32,32,32,32,32,32,32,32,32,32,32,32,125,10,32,32,32,32,32,32,32,32,32,32,32,32,115,99,97,108,101,32,60,49,46,48,48,44,32,49,46,48,48,44,32,49,46,53,48,62,32,42,32,50,46,53,48,10,32,32,32,32,32,32,32,32,32,32,32,32,116,114,97,110,115,108,97,116,101,32,60,48,46,48,48,44,32,48,46,48,48,44,32,48,46,48,48,62,10,32,32,32,32,32,32,32,32,125,10,32,32,32,32,32,32,32,32,102,105,110,105,115,104,32,123,32,97,109,98,105,101,110,116,32,49,32,100,105,102,102,117,115,101,32,48,32,125,10,32,32,32,32,125,10,32,32,32,32,115,99,97,108,101,32,49,48,48,48,48,10,125,10,47,47,109,105,115,116,10,102,111,103,32,123,32,102,111,103,95,116,121,112,101,32,50,10,32,32,32,32,100,105,115,116,97,110,99,101,32,53,48,10,32,32,32,32,99,111,108,111,114,32,114,103,98,32,60,49,46,48,48,44,32,49,46,48,48,44,32,49,46,48,48,62,32,42,32,48,46,55,53,10,32,32,32,32,102,111,103,95,111,102,102,115,101,116,32,48,46,49,48,10,32,32,32,32,102,111,103,95,97,108,116,32,49,46,53,48,10,32,32,32,32,116,117,114,98,117,108,101,110,99,101,32,49,46,55,53,10,125,10,47,47,103,110,100,10,112,108,97,110,101,32,123,32,60,48,46,48,48,44,32,49,46,48,48,44,32,48,46,48,48,62,44,32,48,10,32,32,32,32,116,101,120,116,117,114,101,32,123,10,32,32,32,32,32,32,32,32,112,105,103,109,101,110,116,123,32,99,111,108,111,114,32,114,103,98,32,60,48,46,50,53,44,32,48,46,52,53,44,32,48,46,48,48,62,32,125,10,32,32,32,32,32,32,32,32,110,111,114,109,97,108,32,123,32,98,117,109,112,115,32,48,46,55,53,32,115,99,97,108,101,32,48,46,48,49,32,125,10,32,32,32,32,32,32,32,32,102,105,110,105,115,104,32,123,32,112,104,111,110,103,32,48,46,49,48,32,125,10,32,32,32,32,125,10,125,10,0,0,0,47,121,108,111,114,98,114,51,47,50,0,0,0,0,0,0,99,104,111,99,111,108,97,116,101,0,0,0,0,0,0,0,47,121,108,111,114,98,114,51,47,49,0,0,0,0,0,0,108,65,114,114,0,0,0,0,47,121,108,103,110,98,117,57,47,57,0,0,0,0,0,0,37,46,53,103,32,37,46,53,103,32,37,46,53,103,32,37,115,99,111,108,111,114,10,0,101,114,114,111,114,32,105,110,32,99,111,108,120,108,97,116,101,40,41,10,0,0,0,0,47,121,108,103,110,98,117,57,47,56,0,0,0,0,0,0,47,121,108,103,110,98,117,57,47,55,0,0,0,0,0,0,99,111,112,112,101,114,0,0,47,121,108,103,110,98,117,57,47,54,0,0,0,0,0,0,47,98,117,103,110,57,47,51,0,0,0,0,0,0,0,0,65,118,97,110,116,71,97,114,100,101,45,68,101,109,105,0,60,97,114,101,97,32,115,104,97,112,101,61,34,112,111,108,121,34,0,0,0,0,0,0,47,121,108,103,110,98,117,57,47,53,0,0,0,0,0,0,98,111,108,100,0,0,0,0,47,121,108,103,110,98,117,57,47,52,0,0,0,0,0,0,40,108,105,98,41,58,112,115,0,0,0,0,0,0,0,0,47,121,108,103,110,98,117,57,47,51,0,0,0,0,0,0,97,118,103,95,100,105,115,116,0,0,0,0,0,0,0,0,47,121,108,103,110,98,117,57,47,50,0,0,0,0,0,0,60,115,116,111,112,32,111,102,102,115,101,116,61,34,48,34,32,115,116,121,108,101,61,34,115,116,111,112,45,99,111,108,111,114,58,0,0,0,0,0,47,115,101,116,108,97,121,101,114,32,123,47,109,97,120,108,97,121,101,114,32,101,120,99,104,32,100,101,102,32,47,99,117,114,108,97,121,101,114,32,101,120,99,104,32,100,101,102,0,0,0,0,0,0,0,0,99,97,109,101,114,97,32,123,32,108,111,99,97,116,105,111,110,32,60,37,46,51,102,32,44,32,37,46,51,102,32,44,32,37,46,51,102,62,10,32,32,32,32,32,32,32,32,32,108,111,111,107,95,97,116,32,32,60,37,46,51,102,32,44,32,37,46,51,102,32,44,32,37,46,51,102,62,10,32,32,32,32,32,32,32,32,32,114,105,103,104,116,32,120,32,42,32,105,109,97,103,101,95,119,105,100,116,104,32,47,32,105,109,97,103,101,95,104,101,105,103,104,116,10,32,32,32,32,32,32,32,32,32,97,110,103,108,101,32,37,46,51,102,10,125,10,0,0,0,0,0,0,37,115,32,37,46,51,102,10,0,0,0,0,0,0,0,0,47,121,108,103,110,98,117,57,47,49,0,0,0,0,0,0,47,121,108,103,110,98,117,56,47,56,0,0,0,0,0,0,107,97,112,112,97,0,0,0,97,115,32,114,101,113,117,105,114,101,100,32,98,121,32,116,104,101,32,45,110,32,102,108,97,103,10,0,0,0,0,0,47,121,108,103,110,98,117,56,47,55,0,0,0,0,0,0,47,121,108,103,110,98,117,56,47,54,0,0,0,0,0,0,47,121,108,103,110,98,117,56,47,53,0,0,0,0,0,0,118,111,114,111,95,109,97,114,103,105,110,0,0,0,0,0,47,121,108,103,110,98,117,56,47,52,0,0,0,0,0,0,47,98,117,103,110,57,47,50,0,0,0,0,0,0,0,0,109,97,120,0,0,0,0,0,40,114,118,32,61,61,32,48,41,32,124,124,32,40,78,68,95,111,114,100,101,114,40,114,118,41,45,78,68,95,111,114,100,101,114,40,118,41,41,42,100,105,114,32,62,32,48,0,47,121,108,103,110,98,117,56,47,51,0,0,0,0,0,0,47,121,108,103,110,98,117,56,47,50,0,0,0,0,0,0,47,121,108,103,110,98,117,56,47,49,0,0,0,0,0,0,47,121,108,103,110,98,117,55,47,55,0,0,0,0,0,0,60,100,101,102,115,62,10,60,114,97,100,105,97,108,71,114,97,100,105,101,110,116,32,105,100,61,34,114,95,37,100,34,32,99,120,61,34,53,48,37,37,34,32,99,121,61,34,53,48,37,37,34,32,114,61,34,55,53,37,37,34,32,102,120,61,34,37,100,37,37,34,32,102,121,61,34,37,100,37,37,34,62,10,0,0,0,0,0,47,108,97,121,101,114,108,101,110,32,108,97,121,101,114,99,111,108,111,114,115,101,113,32,108,101,110,103,116,104,32,100,101,102,0,0,0,0,0,0,47,47,42,42,42,32,98,101,103,105,110,95,103,114,97,112,104,32,37,115,10,0,0,0,47,121,108,103,110,98,117,55,47,54,0,0,0,0,0,0,47,121,108,103,110,98,117,55,47,53,0,0,0,0,0,0,105,117,109,108,0,0,0,0,47,121,108,103,110,98,117,55,47,52,0,0,0,0,0,0,47,121,108,103,110,98,117,55,47,51,0,0,0,0,0,0,101,103,103,0,0,0,0,0,99,111,110,115,116,114,97,105,110,101,100,0,0,0,0,0,114,101,97,100,0,0,0,0,47,121,108,103,110,98,117,55,47,50,0,0,0,0,0,0,67,0,0,0,0,0,0,0,47,121,108,103,110,98,117,55,47,49,0,0,0,0,0,0,47,98,117,103,110,57,47,49,0,0,0,0,0,0,0,0,47,121,108,103,110,98,117,54,47,54,0,0,0,0,0,0,47,121,108,103,110,98,117,54,47,53,0,0,0,0,0,0,47,121,108,103,110,98,117,54,47,52,0,0,0,0,0,0,47,121,108,103,110,98,117,54,47,51,0,0,0,0,0,0,37,99,37,103,44,37,103,0,100,101,102,0,0,0,0,0,47,47,42,42,42,32,101,110,100,95,103,114,97,112,104,10,0,0,0,0,0,0,0,0,47,121,108,103,110,98,117,54,47,50,0,0,0,0,0,0,115,112,108,105,116,115,32,105,110,116,111,32,116,119,111,32,110,97,109,101,32,116,111,107,101,110,115,10,0,0,0,0,47,121,108,103,110,98,117,54,47,49,0,0,0,0,0,0,105,115,105,110,0,0,0,0,47,121,108,103,110,98,117,53,47,53,0,0,0,0,0,0,47,121,108,103,110,98,117,53,47,52,0,0,0,0,0,0,47,121,108,103,110,98,117,53,47,51,0,0,0,0,0,0,47,121,108,103,110,98,117,53,47,50,0,0,0,0,0,0,47,98,117,103,110,56,47,56,0,0,0,0,0,0,0,0,47,121,108,103,110,98,117,53,47,49,0,0,0,0,0,0,47,121,108,103,110,98,117,52,47,52,0,0,0,0,0,0,47,121,108,103,110,98,117,52,47,51,0,0,0,0,0,0,47,121,108,103,110,98,117,52,47,50,0,0,0,0,0,0,32,100,61,34,0,0,0,0,9,93,0,0,0,0,0,0,47,47,42,42,42,32,98,101,103,105,110,95,108,97,121,101,114,58,32,37,115,44,32,37,100,47,37,100,10,0,0,0,47,121,108,103,110,98,117,52,47,49,0,0,0,0,0,0,47,121,108,103,110,98,117,51,47,51,0,0,0,0,0,0,105,113,117,101,115,116,0,0,47,121,108,103,110,98,117,51,47,50,0,0,0,0,0,0,47,121,108,103,110,98,117,51,47,49,0,0,0,0,0,0,47,121,108,103,110,57,47,57,0,0,0,0,0,0,0,0,47,121,108,103,110,57,47,56,0,0,0,0,0,0,0,0,47,98,117,103,110,56,47,55,0,0,0,0,0,0,0,0,47,121,108,103,110,57,47,55,0,0,0,0,0,0,0,0,47,121,108,103,110,57,47,54,0,0,0,0,0,0,0,0,47,121,108,103,110,57,47,53,0,0,0,0,0,0,0,0,47,121,108,103,110,57,47,52,0,0,0,0,0,0,0,0,60,112,97,116,104,0,0,0,9,9,91,46,56,32,46,56,32,46,56,93,0,0,0,0,47,47,42,42,42,32,101,110,100,95,108,97,121,101,114,10,0,0,0,0,0,0,0,0,47,121,108,103,110,57,47,51,0,0,0,0,0,0,0,0,47,121,108,103,110,57,47,50,0,0,0,0,0,0,0,0,105,111,116,97,0,0,0,0,47,121,108,103,110,57,47,49,0,0,0,0,0,0,0,0,47,121,108,103,110,56,47,56,0,0,0,0,0,0,0,0,47,121,108,103,110,56,47,55,0,0,0,0,0,0,0,0,47,121,108,103,110,56,47,54,0,0,0,0,0,0,0,0,47,98,117,103,110,56,47,54,0,0,0,0,0,0,0,0,97,103,114,101,97,100,0,0,47,121,108,103,110,56,47,53,0,0,0,0,0,0,0,0,47,121,108,103,110,56,47,52,0,0,0,0,0,0,0,0,47,121,108,103,110,56,47,51,0,0,0,0,0,0,0,0,47,121,108,103,110,56,47,50,0,0,0,0,0,0,0,0,9,9,91,46,54,32,46,56,32,46,56,93,0,0,0,0,47,47,42,42,42,32,98,101,103,105,110,95,112,97,103,101,10,0,0,0,0,0,0,0,47,121,108,103,110,56,47,49,0,0,0,0,0,0,0,0,47,121,108,103,110,55,47,55,0,0,0,0,0,0,0,0,47,121,108,103,110,55,47,54,0,0,0,0,0,0,0,0,47,121,108,103,110,55,47,53,0,0,0,0,0,0,0,0,47,121,108,103,110,55,47,52,0,0,0,0,0,0,0,0,47,121,108,103,110,55,47,51,0,0,0,0,0,0,0,0,47,98,117,103,110,56,47,53,0,0,0,0,0,0,0,0,47,121,108,103,110,55,47,50,0,0,0,0,0,0,0,0,47,121,108,103,110,55,47,49,0,0,0,0,0,0,0,0,47,121,108,103,110,54,47,54,0,0,0,0,0,0,0,0,47,121,108,103,110,54,47,53,0,0,0,0,0,0,0,0,46,46,47,46,46,47,46,46,47,112,108,117,103,105,110,47,99,111,114,101,47,103,118,114,101,110,100,101,114,95,99,111,114,101,95,115,118,103,46,99,0,0,0,0,0,0,0,0,9,9,91,46,52,32,46,56,32,46,56,93,0,0,0,0,47,47,42,42,42,32,101,110,100,95,112,97,103,101,10,0,47,121,108,103,110,54,47,52,0,0,0,0,0,0,0,0,47,121,108,103,110,54,47,51,0,0,0,0,0,0,0,0,105,110,102,105,110,0,0,0,47,121,108,103,110,54,47,50,0,0,0,0,0,0,0,0,47,121,108,103,110,54,47,49,0,0,0,0,0,0,0,0,47,121,108,103,110,53,47,53,0,0,0,0,0,0,0,0,47,121,108,103,110,53,47,52,0,0,0,0,0,0,0,0,47,98,117,103,110,56,47,52,0,0,0,0,0,0,0,0,47,121,108,103,110,53,47,51,0,0,0,0,0,0,0,0,47,121,108,103,110,53,47,50,0,0,0,0,0,0,0,0,47,121,108,103,110,53,47,49,0,0,0,0,0,0,0,0,47,121,108,103,110,52,47,52,0,0,0,0,0,0,0,0,9,9,91,46,50,32,46,56,32,46,56,93,0,0,0,0,47,47,42,42,42,32,98,101,103,105,110,95,99,108,117,115,116,101,114,10,0,0,0,0,47,121,108,103,110,52,47,51,0,0,0,0,0,0,0,0,47,121,108,103,110,52,47,50,0,0,0,0,0,0,0,0,105,109,97,103,101,0,0,0,47,121,108,103,110,52,47,49,0,0,0,0,0,0,0,0,47,121,108,103,110,51,47,51,0,0,0,0,0,0,0,0,47,121,108,103,110,51,47,50,0,0,0,0,0,0,0,0,47,121,108,103,110,51,47,49,0,0,0,0,0,0,0,0,47,98,117,103,110,56,47,51,0,0,0,0,0,0,0,0,47,115,118,103,47,121,101,108,108,111,119,103,114,101,101,110,0,0,0,0,0,0,0,0,47,115,118,103,47,121,101,108,108,111,119,0,0,0,0,0,47,115,118,103,47,119,104,105,116,101,115,109,111,107,101,0,47,115,118,103,47,119,104,105,116,101,0,0,0,0,0,0,53,44,50,0,0,0,0,0,9,9,91,48,32,48,32,48,93,0,0,0,0,0,0,0,47,47,42,42,42,32,101,110,100,95,99,108,117,115,116,101,114,10,0,0,0,0,0,0,47,115,118,103,47,119,104,101,97,116,0,0,0,0,0,0,47,97,99,99,101,110,116,54,47,51,0,0,0,0,0,0,47,115,118,103,47,118,105,111,108,101,116,0,0,0,0,0,105,103,114,97,118,101,0,0,47,115,118,103,47,116,117,114,113,117,111,105,115,101,0,0,104,97,108,102,0,0,0,0,47,115,118,103,47,116,111,109,97,116,111,0,0,0,0,0,47,115,118,103,47,116,104,105,115,116,108,101,0,0,0,0,47,115,118,103,47,116,101,97,108,0,0,0,0,0,0,0,47,98,117,103,110,56,47,50,0,0,0,0,0,0,0,0,101,109,115,99,114,105,112,116,101,110,58,58,109,101,109,111,114,121,95,118,105,101,119,0,47,115,118,103,47,116,97,110,0,0,0,0,0,0,0,0,47,115,118,103,47,115,116,101,101,108,98,108,117,101,0,0,108,105,109,101,0,0,0,0,32,45,116,97,103,115,32,123,37,100,37,115,37,100,125,0,47,115,118,103,47,115,112,114,105,110,103,103,114,101,101,110,0,0,0,0,0,0,0,0,47,115,118,103,47,115,110,111,119,0,0,0,0,0,0,0,49,44,53,0,0,0,0,0,9,91,9,37,32,108,97,121,101,114,32,99,111,108,111,114,32,115,101,113,117,101,110,99,101,32,45,32,100,97,114,107,101,115,116,32,116,111,32,108,105,103,104,116,101,115,116,0,47,47,42,42,42,32,98,101,103,105,110,95,110,111,100,101,58,32,37,115,10,0,0,0,47,115,118,103,47,115,108,97,116,101,103,114,101,121,0,0,99,104,97,114,116,114,101,117,115,101,0,0,0,0,0,0,112,99,0,0,0,0,0,0,47,115,118,103,47,115,108,97,116,101,103,114,97,121,0,0,105,101,120,99,108,0,0,0,47,115,118,103,47,115,108,97,116,101,98,108,117,101,0,0,115,101,116,104,115,98,0,0,37,115,32,105,115,32,110,111,116,32,97,32,107,110,111,119,110,32,99,111,108,111,114,46,10,0,0,0,0,0,0,0,47,115,118,103,47,115,107,121,98,108,117,101,0,0,0,0,47,115,118,103,47,115,105,108,118,101,114,0,0,0,0,0,99,111,111,108,99,111,112,112,101,114,0,0,0,0,0,0,47,115,118,103,47,115,105,101,110,110,97,0,0,0,0,0,47,98,117,103,110,56,47,49,0,0,0,0,0,0,0,0,60,97,114,101,97,32,115,104,97,112,101,61,34,114,101,99,116,34,0,0,0,0,0,0,47,115,118,103,47,115,101,97,115,104,101,108,108,0,0,0,37,100,32,37,100,32,37,100,32,37,100,32,37,100,32,37,100,32,37,100,32,37,100,32,37,100,32,37,46,49,102,32,37,100,32,37,100,32,37,100,32,37,100,10,0,0,0,0,102,105,108,108,101,100,0,0,47,115,118,103,47,115,101,97,103,114,101,101,110,0,0,0,47,115,118,103,47,115,97,110,100,121,98,114,111,119,110,0,102,97,115,116,0,0,0,0,47,115,118,103,47,115,97,108,109,111,110,0,0,0,0,0,47,108,97,121,101,114,99,111,108,111,114,115,101,113,0,0,47,47,42,42,42,32,101,110,100,95,110,111,100,101,10,0,47,115,118,103,47,115,97,100,100,108,101,98,114,111,119,110,0,0,0,0,0,0,0,0,47,115,118,103,47,114,111,121,97,108,98,108,117,101,0,0,105,99,105,114,99,0,0,0,110,111,100,101,32,112,111,115,105,116,105,111,110,115,32,97,114,101,32,105,103,110,111,114,101,100,32,117,110,108,101,115,115,32,115,116,97,114,116,61,114,97,110,100,111,109,10,0,47,115,118,103,47,114,111,115,121,98,114,111,119,110,0,0,47,115,118,103,47,114,101,100,0,0,0,0,0,0,0,0,47,115,118,103,47,112,117,114,112,108,101,0,0,0,0,0,47,115,118,103,47,112,111,119,100,101,114,98,108,117,101,0,47,98,117,103,110,55,47,55,0,0,0,0,0,0,0,0,115,111,117,114,99,101,0,0,47,115,118,103,47,112,108,117,109,0,0,0,0,0,0,0,47,115,118,103,47,112,105,110,107,0,0,0,0,0,0,0,47,115,118,103,47,112,101,114,117,0,0,0,0,0,0,0,47,115,118,103,47,112,101,97,99,104,112,117,102,102,0,0,34,32,115,116,114,111,107,101,45,111,112,97,99,105,116,121,61,34,37,102,0,0,0,0,47,115,104,111,119,112,97,103,101,32,123,32,125,32,100,101,102,0,0,0,0,0,0,0,47,47,42,42,42,32,98,101,103,105,110,95,101,100,103,101,10,0,0,0,0,0,0,0,47,115,118,103,47,112,97,112,97,121,97,119,104,105,112,0,105,110,118,97,108,105,100,32,97,112,105,32,105,110,32,99,111,110,102,105,103,58,32,37,115,32,37,115,10,0,0,0,47,115,118,103,47,112,97,108,101,118,105,111,108,101,116,114,101,100,0,0,0,0,0,0,47,115,118,103,47,112,97,108,101,116,117,114,113,117,111,105,115,101,0,0,0,0,0,0,105,97,99,117,116,101,0,0,37,108,102,37,108,102,37,108,102,0,0,0,0,0,0,0,47,115,118,103,47,112,97,108,101,103,114,101,101,110,0,0,112,111,105,110,116,0,0,0,115,116,97,114,116,32,112,111,114,116,58,32,40,37,46,53,103,44,32,37,46,53,103,41,44,32,116,97,110,103,101,110,116,32,97,110,103,108,101,58,32,37,46,53,103,44,32,37,115,10,0,0,0,0,0,0,37,37,37,37,66,111,117,110,100,105,110,103,66,111,120,58,32,37,100,32,37,100,32,37,100,32,37,100,0,0,0,0,47,115,118,103,47,112,97,108,101,103,111,108,100,101,110,114,111,100,0,0,0,0,0,0,47,115,118,103,47,111,114,99,104,105,100,0,0,0,0,0,47,98,117,103,110,55,47,54,0,0,0,0,0,0,0,0,47,115,118,103,47,111,114,97,110,103,101,114,101,100,0,0,47,115,118,103,47,111,114,97,110,103,101,0,0,0,0,0,47,115,118,103,47,111,108,105,118,101,100,114,97,98,0,0,47,98,117,103,110,55,47,53,0,0,0,0,0,0,0,0,100,101,109,105,0,0,0,0,47,115,118,103,47,111,108,105,118,101,0,0,0,0,0,0,34,32,115,116,114,111,107,101,45,100,97,115,104,97,114,114,97,121,61,34,37,115,0,0,47,101,110,100,112,97,103,101,32,123,32,115,104,111,119,112,97,103,101,32,125,32,98,105,110,100,32,100,101,102,0,0,47,47,42,42,42,32,101,110,100,95,101,100,103,101,10,0,47,115,118,103,47,111,108,100,108,97,99,101,0,0,0,0,47,115,118,103,47,110,97,118,121,0,0,0,0,0,0,0,115,121,110,116,97,120,32,101,114,114,111,114,32,45,32,98,97,100,108,121,32,102,111,114,109,101,100,32,110,117,109,98,101,114,32,39,37,115,39,32,105,110,32,108,105,110,101,32,37,100,32,111,102,32,37,115,10,0,0,0,0,0,0,0,118,101,99,116,111,114,0,0,104,101,108,108,105,112,0,0,47,115,118,103,47,110,97,118,97,106,111,119,104,105,116,101,0,0,0,0,0,0,0,0,47,115,118,103,47,109,111,99,99,97,115,105,110,0,0,0,47,115,118,103,47,109,105,115,116,121,114,111,115,101,0,0,47,115,118,103,47,109,105,110,116,99,114,101,97,109,0,0,95,65,71,95,115,116,114,100,97,116,97,0,0,0,0,0,47,115,118,103,47,109,105,100,110,105,103,104,116,98,108,117,101,0,0,0,0,0,0,0,47,115,118,103,47,109,101,100,105,117,109,118,105,111,108,101,116,114,101,100,0,0,0,0,47,115,118,103,47,109,101,100,105,117,109,116,117,114,113,117,111,105,115,101,0,0,0,0,47,115,118,103,47,109,101,100,105,117,109,115,112,114,105,110,103,103,114,101,101,110,0,0,34,32,115,116,114,111,107,101,45,119,105,100,116,104,61,34,37,103,0,0,0,0,0,0,9,115,101,116,109,97,116,114,105,120,0,0,0,0,0,0,32,32,32,32,110,111,95,115,104,97,100,111,119,10,0,0,47,115,118,103,47,109,101,100,105,117,109,115,108,97,116,101,98,108,117,101,0,0,0,0,47,115,118,103,47,109,101,100,105,117,109,115,101,97,103,114,101,101,110,0,0,0,0,0,104,101,97,114,116,115,0,0,47,115,118,103,47,109,101,100,105,117,109,112,117,114,112,108,101,0,0,0,0,0,0,0,47,115,118,103,47,109,101,100,105,117,109,111,114,99,104,105,100,0,0,0,0,0,0,0,99,114,111,119,0,0,0,0,47,115,118,103,47,109,101,100,105,117,109,98,108,117,101,0,47,115,118,103,47,109,101,100,105,117,109,97,113,117,97,109])
.concat([97,114,105,110,101,0,0,0,47,98,117,103,110,55,47,52,0,0,0,0,0,0,0,0,47,115,118,103,47,109,97,114,111,111,110,0,0,0,0,0,47,115,118,103,47,109,97,103,101,110,116,97,0,0,0,0,47,115,118,103,47,108,105,110,101,110,0,0,0,0,0,0,47,115,118,103,47,108,105,109,101,103,114,101,101,110,0,0,34,32,115,116,114,111,107,101,61,34,0,0,0,0,0,0,9,48,32,48,32,49,32,48,32,51,54,48,32,97,114,99,0,0,0,0,0,0,0,0,116,101,120,116,32,123,10,32,32,32,32,116,116,102,32,34,37,115,34,44,10,32,32,32,32,34,37,115,34,44,32,37,46,51,102,44,32,37,46,51,102,10,32,32,32,32,37,115,32,32,32,32,37,115,32,32,32,32,37,115,32,32,32,32,37,115,32,32,32,32,37,115,125,10,0,0,0,0,0,0,47,115,118,103,47,108,105,109,101,0,0,0,0,0,0,0,47,115,118,103,47,108,105,103,104,116,121,101,108,108,111,119,0,0,0,0,0,0,0,0,104,97,114,114,0,0,0,0,47,115,118,103,47,108,105,103,104,116,115,116,101,101,108,98,108,117,101,0,0,0,0,0,47,115,118,103,47,108,105,103,104,116,115,108,97,116,101,103,114,101,121,0,0,0,0,0,47,115,118,103,47,108,105,103,104,116,115,108,97,116,101,103,114,97,121,0,0,0,0,0,47,115,118,103,47,108,105,103,104,116,115,107,121,98,108,117,101,0,0,0,0,0,0,0,47,98,117,103,110,55,47,51,0,0,0,0,0,0,0,0,97,103,99,108,111,115,101,0,47,115,118,103,47,108,105,103,104,116,115,101,97,103,114,101,101,110,0,0,0,0,0,0,47,115,118,103,47,108,105,103,104,116,115,97,108,109,111,110,0,0,0,0,0,0,0,0,47,115,118,103,47,108,105,103,104,116,112,105,110,107,0,0,47,115,118,103,47,108,105,103,104,116,103,114,101,121,0,0,9,114,120,32,114,121,32,115,99,97,108,101,0,0,0,0,115,99,97,108,101,32,37,46,51,102,10,0,0,0,0,0,47,115,118,103,47,108,105,103,104,116,103,114,101,101,110,0,47,115,118,103,47,108,105,103,104,116,103,114,97,121,0,0,104,65,114,114,0,0,0,0,47,115,118,103,47,108,105,103,104,116,103,111,108,100,101,110,114,111,100,121,101,108,108,111,119,0,0,0,0,0,0,0,47,115,118,103,47,108,105,103,104,116,99,121,97,110,0,0,47,115,118,103,47,108,105,103,104,116,99,111,114,97,108,0,98,111,111,108,0,0,0,0,47,115,118,103,47,108,105,103,104,116,98,108,117,101,0,0,47,98,117,103,110,55,47,50,0,0,0,0,0,0,0,0,47,115,118,103,47,108,101,109,111,110,99,104,105,102,102,111,110,0,0,0,0,0,0,0,47,115,118,103,47,108,97,119,110,103,114,101,101,110,0,0,47,115,118,103,47,108,97,118,101,110,100,101,114,98,108,117,115,104,0,0,0,0,0,0,47,115,118,103,47,108,97,118,101,110,100,101,114,0,0,0,34,32,102,105,108,108,45,111,112,97,99,105,116,121,61,34,37,102,0,0,0,0,0,0,9,120,32,121,32,116,114,97,110,115,108,97,116,101,0,0,47,47,42,42,42,32,116,101,120,116,112,97,114,97,58,32,37,115,44,32,102,111,110,116,115,105,122,101,32,61,32,37,46,51,102,44,32,102,111,110,116,110,97,109,101,32,61,32,37,115,10,0,0,0,0,0,47,115,118,103,47,107,104,97,107,105,0,0,0,0,0,0,47,115,118,103,47,105,118,111,114,121,0,0,0,0,0,0,103,116,0,0,0,0,0,0,47,115,118,103,47,105,110,100,105,103,111,0,0,0,0,0,47,115,118,103,47,105,110,100,105,97,110,114,101,100,0,0,47,115,118,103,47,104,111,116,112,105,110,107,0,0,0,0,47,115,118,103,47,104,111,110,101,121,100,101,119,0,0,0,47,98,117,103,110,55,47,49,0,0,0,0,0,0,0,0,47,115,118,103,47,103,114,101,121,0,0,0,0,0,0,0,47,115,118,103,47,103,114,101,101,110,121,101,108,108,111,119,0,0,0,0,0,0,0,0,47,115,118,103,47,103,114,101,101,110,0,0,0,0,0,0,38,108,116,59,0,0,0,0,47,115,118,103,47,103,114,97,121,0,0,0,0,0,0,0,117,114,108,40,35,114,95,37,100,41,0,0,0,0,0,0,9,110,101,119,112,97,116,104,0,0,0,0,0,0,0,0,115,112,104,101,114,101,32,123,60,37,57,46,51,102,44,32,37,57,46,51,102,44,32,37,57,46,51,102,62,44,32,49,46,48,10,32,32,32,32,37,115,32,32,32,32,37,115,32,32,32,32,37,115,32,32,32,32,37,115,125,10,0,0,0,47,115,118,103,47,103,111,108,100,101,110,114,111,100,0,0,47,115,118,103,47,103,111,108,100,0,0,0,0,0,0,0,103,101,0,0,0,0,0,0,47,115,118,103,47,103,104,111,115,116,119,104,105,116,101,0,47,115,118,103,47,103,97,105,110,115,98,111,114,111,0,0,116,107,58,116,107,0,0,0,47,115,118,103,47,102,117,99,104,115,105,97,0,0,0,0,47,115,118,103,47,102,111,114,101,115,116,103,114,101,101,110,0,0,0,0,0,0,0,0,47,98,117,103,110,54,47,54,0,0,0,0,0,0,0,0,47,115,118,103,47,102,108,111,114,97,108,119,104,105,116,101,0,0,0,0,0,0,0,0,47,115,118,103,47,102,105,114,101,98,114,105,99,107,0,0,47,115,118,103,47,100,111,100,103,101,114,98,108,117,101,0,47,98,117,103,110,54,47,53,0,0,0,0,0,0,0,0,47,115,118,103,47,100,105,109,103,114,101,121,0,0,0,0,117,114,108,40,35,108,95,37,100,41,0,0,0,0,0,0,9,109,97,116,114,105,120,32,99,117,114,114,101,110,116,109,97,116,114,105,120,0,0,0,116,111,114,117,115,32,123,32,37,46,51,102,44,32,37,46,51,102,10,32,32,32,32,37,115,32,32,32,32,37,115,32,32,32,32,37,115,32,32,32,32,37,115,125,10,0,0,0,47,115,118,103,47,100,105,109,103,114,97,121,0,0,0,0,47,97,99,99,101,110,116,54,47,50,0,0,0,0,0,0,47,115,118,103,47,100,101,101,112,115,107,121,98,108,117,101,0,0,0,0,0,0,0,0,103,97,109,109,97,0,0,0,47,115,118,103,47,100,101,101,112,112,105,110,107,0,0,0,47,115,118,103,47,100,97,114,107,118,105,111,108,101,116,0,47,115,118,103,47,100,97,114,107,116,117,114,113,117,111,105,115,101,0,0,0,0,0,0,47,115,118,103,47,100,97,114,107,115,108,97,116,101,103,114,101,121,0,0,0,0,0,0,69,114,114,111,114,0,0,0,101,109,115,99,114,105,112,116,101,110,58,58,118,97,108,0,47,115,118,103,47,100,97,114,107,115,108,97,116,101,103,114,97,121,0,0,0,0,0,0,47,115,118,103,47,100,97,114,107,115,108,97,116,101,98,108,117,101,0,0,0,0,0,0,103,114,101,101,110,0,0,0,47,115,118,103,47,100,97,114,107,115,101,97,103,114,101,101,110,0,0,0,0,0,0,0,47,115,118,103,47,100,97,114,107,115,97,108,109,111,110,0,32,102,105,108,108,61,34,0,9,47,120,32,101,120,99,104,32,100,101,102,0,0,0,0,47,47,42,42,42,32,101,108,108,105,112,115,101,10,0,0,47,115,118,103,47,100,97,114,107,114,101,100,0,0,0,0,99,97,100,101,116,98,108,117,101,0,0,0,0,0,0,0,112,120,0,0,0,0,0,0,47,115,118,103,47,100,97,114,107,111,114,99,104,105,100,0,102,114,97,115,108,0,0,0,47,115,118,103,47,100,97,114,107,111,114,97,110,103,101,0,99,111,108,111,114,32,37,115,0,0,0,0,0,0,0,0,47,115,118,103,47,100,97,114,107,111,108,105,118,101,103,114,101,101,110,0,0,0,0,0,47,115,118,103,47,100,97,114,107,109,97,103,101,110,116,97,0,0,0,0,0,0,0,0,99,108,101,97,114,0,0,0,47,115,118,103,47,100,97,114,107,107,104,97,107,105,0,0,47,98,117,103,110,54,47,52,0,0,0,0,0,0,0,0,37,115,37,115,32,105,115,32,110,111,116,32,97,32,116,114,111,102,102,32,102,111,110,116,10,0,0,0,0,0,0,0,60,97,114,101,97,32,115,104,97,112,101,61,34,99,105,114,99,108,101,34,0,0,0,0,47,115,118,103,47,100,97,114,107,103,114,101,121,0,0,0,110,32,62,61,32,52,0,0,83,32,0,0,0,0,0,0,47,115,118,103,47,100,97,114,107,103,114,101,101,110,0,0,47,115,118,103,47,100,97,114,107,103,114,97,121,0,0,0,115,118,103,58,115,118,103,0,121,101,115,0,0,0,0,0,47,115,118,103,47,100,97,114,107,103,111,108,100,101,110,114,111,100,0,0,0,0,0,0,34,47,62,10,0,0,0,0,9,47,121,32,101,120,99,104,32,100,101,102,0,0,0,0,10,32,32,32,32,37,115,32,32,32,32,37,115,32,32,32,32,37,115,32,32,32,32,37,115,125,10,0,0,0,0,0,47,115,118,103,47,100,97,114,107,99,121,97,110,0,0,0,47,115,118,103,47,100,97,114,107,98,108,117,101,0,0,0,102,114,97,99,51,52,0,0,47,115,118,103,47,99,121,97,110,0,0,0,0,0,0,0,47,115,118,103,47,99,114,105,109,115,111,110,0,0,0,0,47,115,118,103,47,99,111,114,110,115,105,108,107,0,0,0,110,101,119,46,103,118,0,0,100,101,114,105,118,101,100,0,47,115,118,103,47,99,111,114,110,102,108,111,119,101,114,98,108,117,101,0,0,0,0,0,47,98,117,103,110,54,47,51,0,0,0,0,0,0,0,0,109,105,110,0,0,0,0,0,109,99,108,105,109,105,116,0,47,115,118,103,47,99,111,114,97,108,0,0,0,0,0,0,47,115,118,103,47,99,104,111,99,111,108,97,116,101,0,0,47,115,118,103,47,99,104,97,114,116,114,101,117,115,101,0,99,97,110,110,111,116,32,114,101,97,108,108,111,99,32,116,114,105,115,0,0,0,0,0,47,115,118,103,47,99,97,100,101,116,98,108,117,101,0,0,37,103,44,37,103,32,0,0,9,47,114,120,32,101,120,99,104,32,100,101,102,0,0,0,37,115,10,32,32,32,32,37,115,0,0,0,0,0,0,0,47,115,118,103,47,98,117,114,108,121,119,111,111,100,0,0,47,115,118,103,47,98,114,111,119,110,0,0,0,0,0,0,102,114,97,99,49,52,0,0,47,115,118,103,47,98,108,117,101,118,105,111,108,101,116,0,47,115,118,103,47,98,108,117,101,0,0,0,0,0,0,0,99,105,114,99,108,101,0,0,47,115,118,103,47,98,108,97,110,99,104,101,100,97,108,109,111,110,100,0,0,0,0,0,37,100,32,40,37,46,53,103,44,32,37,46,53,103,41,44,32,40,37,46,53,103,44,32,37,46,53,103,41,10,0,0,99,111,117,108,100,110,39,116,32,111,112,101,110,32,101,112,115,102,32,102,105,108,101,32,37,115,10,0,0,0,0,0,92,71,0,0,0,0,0,0,47,115,118,103,47,98,108,97,99,107,0,0,0,0,0,0,47,98,117,103,110,54,47,50,0,0,0,0,0,0,0,0,115,116,114,105,99,116,32,0,47,115,118,103,47,98,105,115,113,117,101,0,0,0,0,0,47,115,118,103,47,98,101,105,103,101,0,0,0,0,0,0,47,115,118,103,47,97,122,117,114,101,0,0,0,0,0,0,47,98,117,103,110,54,47,49,0,0,0,0,0,0,0,0,47,115,118,103,47,97,113,117,97,109,97,114,105,110,101,0,32,112,111,105,110,116,115,61,34,0,0,0,0,0,0,0,9,47,114,121,32,101,120,99,104,32,100,101,102,0,0,0,60,37,57,46,51,102,44,32,37,57,46,51,102,44,32,37,57,46,51,102,62,0,0,0,99,111,108,111,114,115,99,104,101,109,101,0,0,0,0,0,47,115,118,103,47,97,113,117,97,0,0,0,0,0,0,0,112,115,58,112,115,0,0,0,47,115,118,103,47,97,110,116,105,113,117,101,119,104,105,116,101,0,0,0,0,0,0,0,111,117,116,32,111,102,32,100,121,110,97,109,105,99,32,109,101,109,111,114,121,32,105,110,32,97,97,103,95,103,101,116,95,110,101,120,116,95,98,117,102,102,101,114,40,41,0,0,37,46,48,76,102,0,0,0,102,114,97,99,49,50,0,0,47,115,118,103,47,97,108,105,99,101,98,108,117,101,0,0,47,115,112,101,99,116,114,97,108,57,47,57,0,0,0,0,47,115,112,101,99,116,114,97,108,57,47,56,0,0,0,0,47,115,112,101,99,116,114,97,108,57,47,55,0,0,0,0,109,101,109,111,114,121,32,97,108,108,111,99,97,116,105,111,110,32,102,97,105,108,117,114,101,0,0,0,0,0,0,0,47,115,112,101,99,116,114,97,108,57,47,54,0,0,0,0,47,115,112,101,99,116,114,97,108,57,47,53,0,0,0,0,47,115,112,101,99,116,114,97,108,57,47,52,0,0,0,0,105,111,115,116,114,101,97,109,0,0,0,0,0,0,0,0,47,115,112,101,99,116,114,97,108,57,47,51,0,0,0,0,60,112,111,108,121,108,105,110,101,0,0,0,0,0,0,0,47,101,108,108,105,112,115,101,95,112,97,116,104,32,123,0,112,111,108,121,103,111,110,32,123,32,37,100,44,10,0,0,47,115,112,101,99,116,114,97,108,57,47,50,0,0,0,0,47,115,112,101,99,116,114,97,108,57,47,49,0,0,0,0,102,111,114,97,108,108,0,0,47,115,112,101,99,116,114,97,108,56,47,56,0,0,0,0,47,115,112,101,99,116,114,97,108,56,47,55,0,0,0,0,47,115,112,101,99,116,114,97,108,56,47,54,0,0,0,0,47,115,112,101,99,116,114,97,108,56,47,53,0,0,0,0,47,98,117,103,110,53,47,53,0,0,0,0,0,0,0,0,47,115,112,101,99,116,114,97,108,56,47,52,0,0,0,0,47,115,112,101,99,116,114,97,108,56,47,51,0,0,0,0,47,115,112,101,99,116,114,97,108,56,47,50,0,0,0,0,47,115,112,101,99,116,114,97,108,56,47,49,0,0,0,0,9,9,99,108,111,115,101,112,97,116,104,0,0,0,0,0,32,32,32,32,116,111,108,101,114,97,110,99,101,32,48,46,49,10,32,32,32,32,37,115,32,32,32,32,37,115,32,32,32,32,37,115,32,32,32,32,37,115,125,10,0,0,0,0,47,115,112,101,99,116,114,97,108,55,47,55,0,0,0,0,47,115,112,101,99,116,114,97,108,55,47,54,0,0,0,0,102,110,111,102,0,0,0,0,47,115,112,101,99,116,114,97,108,55,47,53,0,0,0,0,47,115,112,101,99,116,114,97,108,55,47,52,0,0,0,0,47,115,112,101,99,116,114,97,108,55,47,51,0,0,0,0,47,115,112,101,99,116,114,97,108,55,47,50,0,0,0,0,47,98,117,103,110,53,47,52,0,0,0,0,0,0,0,0,97,103,111,112,101,110,0,0,47,115,112,101,99,116,114,97,108,55,47,49,0,0,0,0,47,115,112,101,99,116,114,97,108,54,47,54,0,0,0,0,47,115,112,101,99,116,114,97,108,54,47,53,0,0,0,0,112,105,99,58,112,105,99,0,47,115,112,101,99,116,114,97,108,54,47,52,0,0,0,0,60,33,45,45,32,0,0,0,9,9,112,111,112,32,110,101,103,32,48,32,114,108,105,110,101,116,111,0,0,0,0,0,47,47,42,42,42,32,112,111,108,121,103,111,110,10,0,0,47,115,112,101,99,116,114,97,108,54,47,51,0,0,0,0,47,115,112,101,99,116,114,97,108,54,47,50,0,0,0,0,101,120,105,115,116,0,0,0,47,115,112,101,99,116,114,97,108,54,47,49,0,0,0,0,105,115,109,97,112,58,109,97,112,0,0,0,0,0,0,0,47,115,112,101,99,116,114,97,108,53,47,53,0,0,0,0,47,115,112,101,99,116,114,97,108,53,47,52,0,0,0,0,47,115,112,101,99,116,114,97,108,53,47,51,0,0,0,0,47,98,117,103,110,53,47,51,0,0,0,0,0,0,0,0,47,115,112,101,99,116,114,97,108,53,47,50,0,0,0,0,102,105,103,58,102,105,103,0,47,115,112,101,99,116,114,97,108,53,47,49,0,0,0,0,47,115,112,101,99,116,114,97,108,52,47,52,0,0,0,0,47,115,112,101,99,116,114,97,108,52,47,51,0,0,0,0,121,101,108,108,111,119,103,114,101,101,110,0,0,0,0,0,9,9,48,32,101,120,99,104,32,114,108,105,110,101,116,111,0,0,0,0,0,0,0,0,32,32,32,32,32,32,32,32,116,111,108,101,114,97,110,99,101,32,48,46,48,49,10,32,32,32,32,37,115,32,32,32,32,37,115,32,32,32,32,37,115,32,32,32,32,37,115,125,10,0,0,0,0,0,0,0,47,115,112,101,99,116,114,97,108,52,47,50,0,0,0,0,37,100,32,37,100,32,0,0,47,115,112,101,99,116,114,97,108,52,47,49,0,0,0,0,101,117,114,111,0,0,0,0,47,115,112,101,99,116,114,97,108,51,47,51,0,0,0,0,47,115,112,101,99,116,114,97,108,51,47,50,0,0,0,0,47,115,112,101,99,116,114,97,108,51,47,49,0,0,0,0,47,115,112,101,99,116,114,97,108,49,49,47,57,0,0,0,47,98,117,103,110,53,47,50,0,0,0,0,0,0,0,0,47,115,112,101,99,116,114,97,108,49,49,47,56,0,0,0,47,115,112,101,99,116,114,97,108,49,49,47,55,0,0,0,47,115,112,101,99,116,114,97,108,49,49,47,54,0,0,0,47,115,112,101,99,116,114,97,108,49,49,47,53,0,0,0,103,105,102,58,115,118,103,0,9,9,101,120,99,104,32,48,32,114,108,105,110,101,116,111,0,0,0,0,0,0,0,0,98,95,115,112,108,105,110,101,0,0,0,0,0,0,0,0,47,115,112,101,99,116,114,97,108,49,49,47,52,0,0,0,102,100,112,0,0,0,0,0,47,115,112,101,99,116,114,97,108,49,49,47,51,0,0,0,101,117,109,108,0,0,0,0,47,115,112,101,99,116,114,97,108,49,49,47,50,0,0,0,47,115,112,101,99,116,114,97,108,49,49,47,49,49,0,0,47,115,112,101,99,116,114,97,108,49,49,47,49,48,0,0,115,104,97,112,101,0,0,0,105,110,115,101,116,0,0,0,47,115,112,101,99,116,114,97,108,49,49,47,49,0,0,0,47,98,117,103,110,53,47,49,0,0,0,0,0,0,0,0,115,112,101,99,105,102,105,101,100,32,114,111,111,116,32,110,111,100,101,32,34,37,115,34,32,119,97,115,32,110,111,116,32,102,111,117,110,100,46,0,47,115,112,101,99,116,114,97,108,49,48,47,57,0,0,0,83,112,97,114,115,101,77,97,116,114,105,120,95,105,115,95,115,121,109,109,101,116,114,105,99,40,65,44,32,70,65,76,83,69,41,0,0,0,0,0,47,115,112,101,99,116,114,97,108,49,48,47,56,0,0,0,47,115,112,101,99,116,114,97,108,49,48,47,55,0,0,0,47,115,112,101,99,116,114,97,108,49,48,47,54,0,0,0,119,104,105,116,101,115,109,111,107,101,0,0,0,0,0,0,9,9,50,32,99,111,112,121,0,0,0,0,0,0,0,0,47,47,42,42,42,32,98,101,122,105,101,114,10,0,0,0,47,115,112,101,99,116,114,97,108,49,48,47,53,0,0,0,47,97,99,99,101,110,116,54,47,49,0,0,0,0,0,0,47,115,112,101,99,116,114,97,108,49,48,47,52,0,0,0,101,116,104,0,0,0,0,0,47,115,112,101,99,116,114,97,108,49,48,47,51,0,0,0,47,115,112,101,99,116,114,97,108,49,48,47,50,0,0,0,47,115,112,101,99,116,114,97,108,49,48,47,49,48,0,0,115,102,100,112,32,111,110,108,121,32,115,117,112,112,111,114,116,115,32,115,116,97,114,116,61,114,97,110,100,111,109,10,0,0,0,0,0,0,0,0,47,115,112,101,99,116,114,97,108,49,48,47,49,0,0,0,47,98,117,103,110,52,47,52,0,0,0,0,0,0,0,0,115,116,100,58,58,119,115,116,114,105,110,103,0,0,0,0,47,115,101,116,51,57,47,57,0,0,0,0,0,0,0,0,47,115,101,116,51,57,47,56,0,0,0,0,0,0,0,0,103,114,97,121,0,0,0,0,46,46,47,46,46,47,46,46,47,112,108,117,103,105,110,47,99,111,114,101,47,103,118,114,101,110,100,101,114,95,99,111,114,101,95,116,107,46,99,0,47,115,101,116,51,57,47,55,0,0,0,0,0,0,0,0,119,105,100,116,104,32,62,32,48,0,0,0,0,0,0,0,47,115,101,116,51,57,47,54,0,0,0,0,0,0,0,0,9,9,109,111,118,101,116,111,0,0,0,0,0,0,0,0,112,105,103,109,101,110,116,32,123,32,99,111,108,111,114,32,37,115,32,125,10,0,0,0,47,115,101,116,51,57,47,53,0,0,0,0,0,0,0,0,98,117,114,108,121,119,111,111,100,0,0,0,0,0,0,0,47,115,101,116,51,57,47,52,0,0,0,0,0,0,0,0,101,116,97,0,0,0,0,0,47,115,101,116,51,57,47,51,0,0,0,0,0,0,0,0,47,115,101,116,51,57,47,50,0,0,0,0,0,0,0,0,47,115,101,116,51,57,47,49,0,0,0,0,0,0,0,0,47,115,101,116,51,56,47,56,0,0,0,0,0,0,0,0,47,98,117,103,110,52,47,51,0,0,0,0,0,0,0,0,37,48,51,111,0,0,0,0,108,101,110,0,0,0,0,0,114,101,99,116,97,110,103,108,101,32,40,37,100,44,37,100,41,32,40,37,100,44,37,100,41,32,37,115,32,37,115,10,0,0,0,0,0,0,0,0,47,115,101,116,51,56,47,55,0,0,0,0,0,0,0,0,46,46,47,46,46,47,46,46,47,112,108,117,103,105,110,47,99,111,114,101,47,103,118,114,101,110,100,101,114,95,99,111,114,101,95,102,105,103,46,99,0,0,0,0,0,0,0,0,115,101,116,108,105,110,101,119,105,100,116,104,40,37,46,51,102,41,0,0,0,0,0,0,47,115,101,116,51,56,47,54,0,0,0,0,0,0,0,0,106,112,103,58,118,114,109,108,0,0,0,0,0,0,0,0,47,115,101,116,51,56,47,53,0,0,0,0,0,0,0,0,120,33,61,78,85,76,76,0,116,114,117,101,0,0,0,0,47,115,101,116,51,56,47,52,0,0,0,0,0,0,0,0,119,104,101,97,116,0,0,0,9,9,52,32,50,32,114,111,108,108,0,0,0,0,0,0,47,115,101,116,51,56,47,51,0,0,0,0,0,0,0,0,47,115,101,116,51,56,47,50,0,0,0,0,0,0,0,0,101,113,117,105,118,0,0,0,114,97,110,100,111,109,0,0,47,115,101,116,51,56,47,49,0,0,0,0,0,0,0,0,47,115,101,116,51,55,47,55,0,0,0,0,0,0,0,0,112,111,108,121,103,111,110,0,47,115,101,116,51,55,47,54,0,0,0,0,0,0,0,0,47,115,101,116,51,55,47,53,0,0,0,0,0,0,0,0,78,111,32,108,111,97,100,105,109,97,103,101,32,112,108,117,103,105,110,32,102,111,114,32,34,37,115,34,10,0,0,0,47,98,117,103,110,52,47,50,0,0,0,0,0,0,0,0,114,97,110,107,0,0,0,0,47,115,101,116,51,55,47,52,0,0,0,0,0,0,0,0,116,97,105,108,112,111,114,116,0,0,0,0,0,0,0,0,47,115,101,116,51,55,47,51,0,0,0,0,0,0,0,0,47,115,101,116,51,55,47,50,0,0,0,0,0,0,0,0,99,97,110,110,111,116,32,109,97,108,108,111,99,32,116,114,105,115,0,0,0,0,0,0,47,115,101,116,51,55,47,49,0,0,0,0,0,0,0,0,118,105,111,108,101,116,0,0,47,98,111,120,112,114,105,109,32,123,9,9,9,9,37,32,120,99,111,114,110,101,114,32,121,99,111,114,110,101,114,32,120,115,105,122,101,32,121,115,105,122,101,0,0,0,0,0,46,46,47,46,46,47,46,46,47,112,108,117,103,105,110,47,99,111,114,101,47,103,118,114,101,110,100,101,114,95,99,111,114,101,95,112,111,118,46,99,0,0,0,0,0,0,0,0,47,115,101,116,51,54,47,54,0,0,0,0,0,0,0,0,37,115,32,105,115,32,122,101,114,111,32,115,105,122,101,100,44,32,111,114,32,111,116,104,101,114,32,114,101,97,100,32,101,114,114,111,114,46,10,0,47,115,101,116,51,54,47,53,0,0,0,0,0,0,0,0,47,115,101,116,51,54,47,52,0,0,0,0,0,0,0,0,47,115,101,116,51,54,47,51,0,0,0,0,0,0,0,0,111,118,97,108,0,0,0,0,47,115,101,116,51,54,47,50,0,0,0,0,0,0,0,0,37,100,32,98,111,120,101,115,58,10,0,0,0,0,0,0,85,84,70,45,56,32,105,110,112,117,116,32,117,115,101,115,32,110,111,110,45,76,97,116,105,110,49,32,99,104,97,114,97,99,116,101,114,115,32,119,104,105,99,104,32,99,97,110,110,111,116,32,98,101,32,104,97,110,100,108,101,100,32,98,121,32,116,104,105,115,32,80,111,115,116,83,99,114,105,112,116,32,100,114,105,118,101,114,10,0,0,0,0,0,0,0,37,108,102,44,37,108,102,44,37,108,102,37,99,0,0,0,47,115,101,116,51,54,47,49,0,0,0,0,0,0,0,0,47,98,117,103,110,52,47,49,0,0,0,0,0,0,0,0,47,115,101,116,51,53,47,53,0,0,0,0,0,0,0,0,100,105,0,0,0,0,0,0,47,115,101,116,51,53,47,52,0,0,0,0,0,0,0,0,47,115,101,116,51,53,47,51,0,0,0,0,0,0,0,0,105,116,97,108,105,99,0,0,47,115,101,116,51,53,47,50,0,0,0,0,0,0,0,0,116,117,114,113,117,111,105,115,101,0,0,0,0,0,0,0,9,103,114,101,115,116,111,114,101,0,0,0,0,0,0,0,111,111,112,115,44,32,105,110,116,101,114,110,97,108,32,101,114,114,111,114,58,32,117,110,104,97,110,100,108,101,100,32,99,111,108,111,114,32,116,121,112,101,61,37,100,32,37,115,10,0,0,0,0,0,0,0,47,115,101,116,51,53,47,49,0,0,0,0,0,0,0,0,47,115,101,116,51,52,47,52,0,0,0,0,0,0,0,0,105,110,112,117,116,32,105,110,32,102,108,101,120,32,115,99,97,110,110,101,114,32,102,97,105,108,101,100,0,0,0,0,109,111,110,101,121,95,103,101,116,32,101,114,114,111,114,0,101,110,115,112,0,0,0,0,47,115,101,116,51,52,47,51,0,0,0,0,0,0,0,0,47,115,101,116,51,52,47,50,0,0,0,0,0,0,0,0,47,115,101,116,51,52,47,49,0,0,0,0,0,0,0,0,47,115,101,116,51,51,47,51,0,0,0,0,0,0,0,0,47,98,117,103,110,51,47,51,0,0,0,0,0,0,0,0,47,115,101,116,51,51,47,50,0,0,0,0,0,0,0,0,47,115,101,116,51,51,47,49,0,0,0,0,0,0,0,0,47,115,101,116,51,49,50,47,57,0,0,0,0,0,0,0,47,115,101,116,51,49,50,47,56,0,0,0,0,0,0,0,100,101,108,97,117,110,97,121,95,116,114,105,97,110,103,117,108,97,116,105,111,110,58,32,37,115,10,0,0,0,0,0,116,111,109,97,116,111,0,0,9,9,125,32,105,102,0,0,114,103,98,60,37,57,46,51,102,44,32,37,57,46,51,102,44,32,37,57,46,51,102,62,32,116,114,97,110,115,109,105,116,32,37,46,51,102,0,0,47,115,101,116,51,49,50,47,55,0,0,0,0,0,0,0,46,46,47,46,46,47,46,46,47,108,105,98,47,110,101,97,116,111,103,101,110,47,99,111,110,115,116,114,97,105,110,116,46,99,0,0,0,0,0,0,47,115,101,116,51,49,50,47,54,0,0,0,0,0,0,0,101,109,115,112,0,0,0,0,47,115,101,116,51,49,50,47,53,0,0,0,0,0,0,0,47,115,101,116,51,49,50,47,52,0,0,0,0,0,0,0,47,115,101,116,51,49,50,47,51,0,0,0,0,0,0,0,47,115,101,116,51,49,50,47,50,0,0,0,0,0,0,0,47,98,117,103,110,51,47,50,0,0,0,0,0,0,0,0,47,115,101,116,51,49,50,47,49,50,0,0,0,0,0,0,47,115,101,116,51,49,50,47,49,49,0,0,0,0,0,0,119,101,105,103,104,116,0,0,47,115,101,116,51,49,50,47,49,48,0,0,0,0,0,0,47,115,101,116,51,49,50,47,49,0,0,0,0,0,0,0,41,10,45,45,62,10,0,0,116,104,105,115,116,108,101,0,9,9,9,116,101,120,116,32,115,116,114,105,110,103,119,105,100,116,104,32,112,111,112,32,119,105,100,116,104,32,101,120,99,104,32,115,117,98,32,116,101,120,116,32,108,101,110,103,116,104,32,100,105,118,32,48,32,116,101,120,116,32,97,115,104,111,119,0,0,0,0,0,66,108,117,101,0,0,0,0,47,115,101,116,51,49,49,47,57,0,0,0,0,0,0,0,47,115,101,116,51,49,49,47,56,0,0,0,0,0,0,0,101,109,112,116,121,0,0,0,47,115,101,116,51,49,49,47,55,0,0,0,0,0,0,0,119,105,100,116,104,0,0,0,47,115,101,116,51,49,49,47,54,0,0,0,0,0,0,0,47,115,101,116,51,49,49,47,53,0,0,0,0,0,0,0,115,97,109,101,116,97,105,108,0,0,0,0,0,0,0,0,47,115,101,116,51,49,49,47,52,0,0,0,0,0,0,0,47,98,117,103,110,51,47,49,0,0,0,0,0,0,0,0,103,118,70,114,101,101,67,111,110,116,101,120,116,0,0,0,47,115,101,116,51,49,49,47,51,0,0,0,0,0,0,0,47,115,101,116,51,49,49,47,50,0,0,0,0,0,0,0,47,115,101,116,51,49,49,47,49,49,0,0,0,0,0,0,110,101,119,114,97,110,107,0,47,115,101,116,51,49,49,47,49,48,0,0,0,0,0,0,32,40,0,0,0,0,0,0,9,9,9,91,93,32,48,32,115,101,116,100,97,115,104,0,71,114,101,101,110,0,0,0,47,115,101,116,51,49,49,47,49,0,0,0,0,0,0,0,114,97,110,107,40,103,44,32,50,44,32,110,115,105,116,101,114,50,40,103,41,41,32,61,61,32,48,0,0,0,0,0,47,115,101,116,51,49,48,47,57,0,0,0,0,0,0,0,47,115,101,116,51,49,48,47,56,0,0,0,0,0,0,0,101,103,114,97,118,101,0,0,76,101,102,116,0,0,0,0,47,115,101,116,51,49,48,47,55,0,0,0,0,0,0,0,47,115,101,116,51,49,48,47,54,0,0,0,0,0,0,0,105,110,115,116,97,108,108,95,105,110,95,114,97,110,107,44,32,108,105,110,101,32,37,100,58,32,37,115,32,37,115,32,114,97,110,107,32,37,100,32,105,32,61,32,37,100,32,97,110,32,61,32,48,10,0,0,78,68,95,111,117,116,40,118,41,46,115,105,122,101,32,61,61,32,50,0,0,0,0,0,47,115,101,116,51,49,48,47,53,0,0,0,0,0,0,0,47,98,114,98,103,57,47,57,0,0,0,0,0,0,0,0,47,115,101,116,51,49,48,47,52,0,0,0,0,0,0,0,47,115,101,116,51,49,48,47,51,0,0,0,0,0,0,0,101,32,33,61,32,78,85,76,76,0,0,0,0,0,0,0,47,115,101,116,51,49,48,47,50,0,0,0,0,0,0,0,47,115,101,116,51,49,48,47,49,48,0,0,0,0,0,0,32,118,101,114,115,105,111,110,32,0,0,0,0,0,0,0,116,97,110,0,0,0,0,0,9,9,119,105,100,116,104,32,48,32,103,116,32,123,0,0,82,101,100,0,0,0,0,0,47,115,101,116,51,49,48,47,49,0,0,0,0,0,0,0,46,46,47,46,46,47,46,46,47,108,105,98,47,100,111,116,103,101,110,47,100,111,116,115,112,108,105,110,101,115,46,99,0,0,0,0,0,0,0,0,47,115,101,116,50,56,47,56,0,0,0,0,0,0,0,0,101,99,105,114,99,0,0,0,47,115,101,116,50,56,47,55,0,0,0,0,0,0,0,0,100,111,116,32,100,111,101,115,32,110,111,116,32,115,117,112,112,111,114,116,32,116,104,101,32,97,115,112,101,99,116,32,97,116,116,114,105,98,117,116,101,32,102,111,114,32,100,105,115,99,111,110,110,101,99,116,101,100,32,103,114,97,112,104,115,32,111,114,32,103,114,97,112,104,115,32,119,105,116,104,32,99,108,117,115,116,101,114,115,10,0,0,0,0,0,0,114,101,98,117,105,108,116,100,95,118,108,105,115,116,115,58,32,114,97,110,107,32,108,101,97,100,32,37,115,32,110,111,116,32,105,110,32,111,114,100,101,114,32,37,100,32,111,102,32,114,97,110,107,32,37,100,10,0,0,0,0,0,0,0,47,115,101,116,50,56,47,54,0,0,0,0,0,0,0,0,47,115,101,116,50,56,47,53,0,0,0,0,0,0,0,0,108,116,97,105,108,0,0,0,47,115,101,116,50,56,47,52,0,0,0,0,0,0,0,0,47,98,114,98,103,57,47,56,0,0,0,0,0,0,0,0,46,46,47,46,46,47,46,46,47,108,105,98,47,100,111,116,103,101,110,47,99,108,117,115,116,101,114,46,99,0,0,0,47,115,101,116,50,56,47,51,0,0,0,0,0,0,0,0,47,115,101,116,50,56,47,50,0,0,0,0,0,0,0,0,47,115,101,116,50,56,47,49,0,0,0,0,0,0,0,0,67,111,117,108,100,32,110,111,116,32,111,112,101,110,32,34,37,115,34,32,102,111,114,32,119,114,105,116,105,110,103,32,58,32,37,115,10,0,0,0,97,99,116,117,97,108,0,0,47,115,101,116,50,55,47,55,0,0,0,0,0,0,0,0,10,60,33,45,45,32,71,101,110,101,114,97,116,101,100,32,98,121,32,0,0,0,0,0,115,116,101,101,108,98,108,117,101,0,0,0,0,0,0,0,9,103,115,97,118,101,0,0,37,115,32,116,114,97,110,115,109,105,116,32,37,46,51,102,0,0,0,0,0,0,0,0,47,115,101,116,50,55,47,54,0,0,0,0,0,0,0,0,47,115,101,116,50,55,47,53,0,0,0,0,0,0,0,0,95,98,108,111,99,107,95,37,100,0,0,0,0,0,0,0,101,97,99,117,116,101,0,0,47,115,101,116,50,55,47,52,0,0,0,0,0,0,0,0,47,115,101,116,50,55,47,51,0,0,0,0,0,0,0,0,47,115,101,116,50,55,47,50,0,0,0,0,0,0,0,0,47,115,101,116,50,55,47,49,0,0,0,0,0,0,0,0,50,46,51,48,46,49,0,0,47,98,114,98,103,57,47,55,0,0,0,0,0,0,0,0,47,115,101,116,50,54,47,54,0,0,0,0,0,0,0,0,46,46,47,46,46,47,46,46,47,108,105,98,47,112,97,116,104,112,108,97,110,47,115,104,111,114,116,101,115,116,46,99,0,0,0,0,0,0,0,0,47,115,101,116,50,54,47,53,0,0,0,0,0,0,0,0,46,46,47,46,46,47,46,46,47,108,105,98,47,112,97,116,104,112,108,97,110,47,114,111,117,116,101,46,99,0,0,0,106,32,61,61,32,48,0,0,47,115,101,116,50,54,47,52,0,0,0,0,0,0,0,0,47,115,101,116,50,54,47,51,0,0,0,0,0,0,0,0,60,72,84,77,76,62,10,0,115,112,114,105,110,103,103,114,101,101,110,0,0,0,0,0,9,47,119,105,100,116,104,32,101,120,99,104,32,100,101,102,0,0,0,0,0,0,0,0,37,115,37,115,0,0,0,0,47,115,101,116,50,54,47,50,0,0,0,0,0,0,0,0,47,97,99,99,101,110,116,53,47,53,0,0,0,0,0,0,47,115,101,116,50,54,47,49,0,0,0,0,0,0,0,0,100,105,118,105,100,101,0,0,47,115,101,116,50,53,47,53,0,0,0,0,0,0,0,0,47,115,101,116,50,53,47,52,0,0,0,0,0,0,0,0,47,115,101,116,50,53,47,51,0,0,0,0,0,0,0,0,47,115,101,116,50,53,47,50,0,0,0,0,0,0,0,0,47,98,114,98,103,57,47,54,0,0,0,0,0,0,0,0,115,116,100,58,58,115,116,114,105,110,103,0,0,0,0,0,47,115,101,116,50,53,47,49,0,0,0,0,0,0,0,0,47,115,101,116,50,52,47,52,0,0,0,0,0,0,0,0,102,117,99,104,115,105,97,0,103,114,97,112,104,32,108,97,98,101,108,0,0,0,0,0,47,115,101,116,50,52,47,51,0,0,0,0,0,0,0,0,47,115,101,116,50,52,47,50,0,0,0,0,0,0,0,0,32,99,111,111,114,100,111,114,105,103,105,110,61,34,48,44,48,34,32,99,111,111,114,100,115,105,122,101,61,34,37,100,44,37,100,34,32,62,0,0,115,110,111,119,0,0,0,0,9,47,116,101,120,116,32,101,120,99,104,32,100,101,102,0,32,32,32,32,116,111,108,101,114,97,110,99,101,32,48,46,48,49,10,32,32,32,32,37,115,32,32,32,32,37,115,32,32,32,32,37,115,32,32,32,32,37,115,125,10,0,0,0,47,115,101,116,50,52,47,49,0,0,0,0,0,0,0,0,98,114,111,119,110,0,0,0,37,108,102,32,37,108,102,32,37,108,102,32,37,108,102,0,47,115,101,116,50,51,47,51,0,0,0,0,0,0,0,0,100,105,97,109,115,0,0,0,47,115,101,116,50,51,47,50,0,0,0,0,0,0,0,0,47,115,101,116,50,51,47,49,0,0,0,0,0,0,0,0,47,115,101,116,49,57,47,57,0,0,0,0,0,0,0,0,47,115,101,116,49,57,47,56,0,0,0,0,0,0,0,0,47,98,114,98,103,57,47,53,0,0,0,0,0,0,0,0,34,37,115,34,32,97,116,32,40,37,46,53,102,44,37,46,53,102,41,59,10,0,0,0,47,115,101,116,49,57,47,55,0,0,0,0,0,0,0,0,35,37,48,50,120,37,48,50,120,37,48,50,120,37,48,50,120,0,0,0,0,0,0,0,47,115,101,116,49,57,47,54,0,0,0,0,0,0,0,0,106,112,101,58,118,114,109,108,0,0,0,0,0,0,0,0,47,115,101,116,49,57,47,53,0,0,0,0,0,0,0,0,110,111,114,109,97,108,0,0,47,115,101,116,49,57,47,52,0,0,0,0,0,0,0,0,32,119,105,100,116,104,58,32,37,100,112,116,59,32,104,101,105,103,104,116,58,32,37,100,112,116,34,0,0,0,0,0,115,108,97,116,101,103,114,101,121,0,0,0,0,0,0,0,47,97,108,105,103,110,101,100,116,101,120,116,32,123,9,9,9,37,32,119,105,100,116,104,32,116,101,120,116,0,0,0,37,115,32,32,32,32,37,115,0,0,0,0,0,0,0,0,77,97,120,46,32,105,116,101,114,97,116,105,111,110,115,32,40,37,100,41,32,114,101,97,99,104,101,100,32,111,110,32,103,114,97,112,104,32,37,115,10,0,0,0,0,0,0,0,47,115,101,116,49,57,47,51,0,0,0,0,0,0,0,0,47,115,101,116,49,57,47,50,0,0,0,0,0,0,0,0,100,101,108,116,97,0,0,0,114,101,103,117,108,97,114,0,47,115,101,116,49,57,47,49,0,0,0,0,0,0,0,0,83,97,116,0,0,0,0,0,47,115,101,116,49,56,47,56,0,0,0,0,0,0,0,0,47,115,101,116,49,56,47,55,0,0,0,0,0,0,0,0,37,102,44,37,102,0,0,0,95,100,103,95,37,100,0,0,47,115,101,116,49,56,47,54,0,0,0,0,0,0,0,0,47,98,114,98,103,57,47,52,0,0,0,0,0,0,0,0,95,110,101,119,95,114,97,110,107,0,0,0,0,0,0,0,47,115,101,116,49,56,47,53,0,0,0,0,0,0,0,0,104,101,97,100,112,111,114,116,0,0,0,0,0,0,0,0,99,108,117,115,116,101,114,32,110,97,109,101,100,32,37,115,32,110,111,116,32,102,111,117,110,100,10,0,0,0,0,0,47,115,101,116,49,56,47,52,0,0,0,0,0,0,0,0,47,115,101,116,49,56,47,51,0,0,0,0,0,0,0,0,116,114,105,97,110,103,117,108,97,116,105,111,110,32,102,97,105,108,101,100,0,0,0,0,70,114,105,0,0,0,0,0,47,115,101,116,49,56,47,50,0,0,0,0,0,0,0,0,32,60,118,58,103,114,111,117,112,32,115,116,121,108,101,61,34,112,111,115,105,116,105,111,110,58,114,101,108,97,116,105,118,101,59,32,0,0,0,0,115,108,97,116,101,103,114,97,121,0,0,0,0,0,0,0,37,32,100,114,97,119,32,116,101,120,116,32,102,105,116,116,101,100,32,116,111,32,105,116,115,32,101,120,112,101,99,116,101,100,32,119,105,100,116,104,0,0,0,0,0,0,0,0,60,37,57,46,51,102,44,32,37,57,46,51,102,44,32,37,57,46,51,102,62,44,32,37,46,51,102,10,0,0,0,0,47,115,101,116,49,56,47,49,0,0,0,0,0,0,0,0,102,97,105,108,101,100,32,116,111,32,111,112,101,110,32,37,115,32,102,111,114,32,114,101,97,100,46,10,0,0,0,0,47,115,101,116,49,55,47,55,0,0,0,0,0,0,0,0,100,101,103,0,0,0,0,0,47,115,101,116,49,55,47,54,0,0,0,0,0,0,0,0,47,115,101,116,49,55,47,53,0,0,0,0,0,0,0,0,125,32,98,105,110,100,32,100,101,102,10,0,0,0,0,0,47,115,101,116,49,55,47,52,0,0,0,0,0,0,0,0,47,115,101,116,49,55,47,51,0,0,0,0,0,0,0,0,47,98,114,98,103,57,47,51,0,0,0,0,0,0,0,0,47,115,101,116,49,55,47,50,0,0,0,0,0,0,0,0,67,58,47,77,105,110,71,87,47,109,115,121,115,47,49,46,48,47,108,111,99,97,108,47,108,105,98,47,103,114,97,112,104,118,105,122,0,0,0,0,47,115,101,116,49,55,47,49,0,0,0,0,0,0,0,0,47,115,101,116,49,54,47,54,0,0,0,0,0,0,0,0,111,98,108,105,113,117,101,0,84,104,117,0,0,0,0,0,47,115,101,116,49,54,47,53,0,0,0,0,0,0,0,0,60,120,109,108,58,110,97,109,101,115,112,97,99,101,32,110,115,61,34,117,114,110,58,115,99,104,101,109,97,115,45,109,105,99,114,111,115,111,102,116,45,99,111,109,58,118,109,108,34,32,112,114,101,102,105,120,61,34,118,34,32,47,62,10,0,0,0,0,0,0,0,0,115,108,97,116,101,98,108,117,101,0,0,0,0,0,0,0,9,115,99,97,108,101,102,111,110,116,32,115,101,116,102,111,110,116,0,0,0,0,0,0,108,105,110,101,97,114,95,115,112,108,105,110,101,0,0,0,47,115,101,116,49,54,47,52,0,0,0,0,0,0,0,0,47,115,101,116,49,54,47,51,0,0,0,0,0,0,0,0,102,97,116,97,108,32,101,114,114,111,114,32,45,32,115,99,97,110,110,101,114,32,105,110,112,117,116,32,98,117,102,102,101,114,32,111,118,101,114,102,108,111,119,0,0,0,0,0,37,76,102,0,0,0,0,0,100,97,114,114,0,0,0,0,47,115,101,116,49,54,47,50,0,0,0,0,0,0,0,0,85,110,115,117,112,112,111,114,116,101,100,32,99,104,97,114,115,101,116,32,34,37,115,34,32,45,32,97,115,115,117,109,105,110,103,32,117,116,102,45,56,10,0,0,0,0,0,0,47,115,101,116,49,54,47,49,0,0,0,0,0,0,0,0,47,115,101,116,49,53,47,53,0,0,0,0,0,0,0,0,47,115,101,116,49,53,47,52,0,0,0,0,0,0,0,0,47,98,114,98,103,57,47,50,0,0,0,0,0,0,0,0,47,115,101,116,49,53,47,51,0,0,0,0,0,0,0,0,47,115,101,116,49,53,47,50,0,0,0,0,0,0,0,0,47,115,101,116,49,53,47,49,0,0,0,0,0,0,0,0,46,46,47,46,46,47,46,46,47,108,105,98,47,99,111,109,109,111,110,47,115,104,97,112,101,115,46,99,0,0,0,0,87,101,100,0,0,0,0,0,47,115,101,116,49,52,47,52,0,0,0,0,0,0,0,0,60,47,83,84,89,76,69,62,10,0,0,0,0,0,0,0,115,107,121,98,108,117,101,0,9,102,105,110,100,102,111,110])
.concat([116,32,101,120,99,104,0,0,115,112,104,101,114,101,95,115,119,101,101,112,32,123,10,32,32,32,32,37,115,10,32,32,32,32,37,100,44,10,0,0,47,115,101,116,49,52,47,51,0,0,0,0,0,0,0,0,47,115,101,116,49,52,47,50,0,0,0,0,0,0,0,0,100,97,103,103,101,114,0,0,47,115,101,116,49,52,47,49,0,0,0,0,0,0,0,0,117,116,102,56,0,0,0,0,115,104,97,112,101,102,105,108,101,32,110,111,116,32,115,101,116,32,111,114,32,110,111,116,32,102,111,117,110,100,32,102,111,114,32,101,112,115,102,32,110,111,100,101,32,37,115,10,0,0,0,0,0,0,0,0,47,115,101,116,49,51,47,51,0,0,0,0,0,0,0,0,47,115,101,116,49,51,47,50,0,0,0,0,0,0,0,0,47,115,101,116,49,51,47,49,0,0,0,0,0,0,0,0,47,98,114,98,103,57,47,49,0,0,0,0,0,0,0,0,47,114,101,100,115,57,47,57,0,0,0,0,0,0,0,0,47,114,101,100,115,57,47,56,0,0,0,0,0,0,0,0,47,114,101,100,115,57,47,55,0,0,0,0,0,0,0,0,84,117,101,0,0,0,0,0,47,114,101,100,115,57,47,54,0,0,0,0,0,0,0,0,118,92,58,42,32,123,32,98,101,104,97,118,105,111,114,58,32,117,114,108,40,35,100,101,102,97,117,108,116,35,86,77,76,41,59,100,105,115,112,108,97,121,58,105,110,108,105,110,101,45,98,108,111,99,107,125,10,0,0,0,0,0,0,0,47,115,101,116,95,102,111,110,116,32,123,0,0,0,0,0,116,114,97,110,115,108,97,116,101,60,37,57,46,51,102,44,32,37,57,46,51,102,44,32,37,57,46,51,102,62,10,0,47,114,101,100,115,57,47,53,0,0,0,0,0,0,0,0,47,114,101,100,115,57,47,52,0,0,0,0,0,0,0,0,100,65,114,114,0,0,0,0,47,114,101,100,115,57,47,51,0,0,0,0,0,0,0,0,98,105,103,53,0,0,0,0,47,114,101,100,115,57,47,50,0,0,0,0,0,0,0,0,105,110,32,108,97,98,101,108,32,111,102,32,110,111,100,101,32,37,115,10,0,0,0,0,47,114,101,100,115,57,47,49,0,0,0,0,0,0,0,0,47,114,101,100,115,56,47,56,0,0,0,0,0,0,0,0,47,98,114,98,103,56,47,56,0,0,0,0,0,0,0,0,103,118,70,114,101,101,76,97,121,111,117,116,0,0,0,0,47,114,101,100,115,56,47,55,0,0,0,0,0,0,0,0,47,114,101,100,115,56,47,54,0,0,0,0,0,0,0,0,47,114,101,100,115,56,47,53,0,0,0,0,0,0,0,0,77,111,110,0,0,0,0,0,47,114,101,100,115,56,47,52,0,0,0,0,0,0,0,0,60,83,84,89,76,69,62,10,0,0,0,0,0,0,0,0,115,105,101,110,110,97,0,0,9,125,32,105,102,0,0,0,114,111,116,97,116,101,32,32,32,60,37,57,46,51,102,44,32,37,57,46,51,102,44,32,37,57,46,51,102,62,10,0,47,114,101,100,115,56,47,51,0,0,0,0,0,0,0,0,47,114,101,100,115,56,47,50,0,0,0,0,0,0,0,0,99,117,114,114,101,110,0,0,47,114,101,100,115,56,47,49,0,0,0,0,0,0,0,0,98,105,103,45,53,0,0,0,47,114,101,100,115,55,47,55,0,0,0,0,0,0,0,0,47,114,101,100,115,55,47,54,0,0,0,0,0,0,0,0,47,114,101,100,115,55,47,53,0,0,0,0,0,0,0,0,47,98,114,98,103,56,47,55,0,0,0,0,0,0,0,0,97,108,112,104,97,0,0,0,47,114,101,100,115,55,47,52,0,0,0,0,0,0,0,0,47,114,101,100,115,55,47,51,0,0,0,0,0,0,0,0,47,114,101,100,115,55,47,50,0,0,0,0,0,0,0,0,47,114,101,100,115,55,47,49,0,0,0,0,0,0,0,0,32,119,105,100,116,104,58,32,37,100,112,116,59,32,104,101,105,103,104,116,58,32,37,100,112,116,34,62,10,0,0,0,115,101,97,115,104,101,108,108,0,0,0,0,0,0,0,0,9,9,103,114,101,115,116,111,114,101,0,0,0,0,0,0,115,99,97,108,101,32,32,32,32,60,37,57,46,51,102,44,32,37,57,46,51,102,44,32,37,57,46,51,102,62,10,0,47,114,101,100,115,54,47,54,0,0,0,0,0,0,0,0,47,114,101,100,115,54,47,53,0,0,0,0,0,0,0,0,99,117,112,0,0,0,0,0,47,114,101,100,115,54,47,52,0,0,0,0,0,0,0,0,73,83,79,45,73,82,45,49,48,48,0,0,0,0,0,0,83,117,110,0,0,0,0,0,47,114,101,100,115,54,47,51,0,0,0,0,0,0,0,0,47,114,101,100,115,54,47,50,0,0,0,0,0,0,0,0,47,114,101,100,115,54,47,49,0,0,0,0,0,0,0,0,47,98,114,98,103,56,47,54,0,0,0,0,0,0,0,0,47,114,101,100,115,53,47,53,0,0,0,0,0,0,0,0,78,111,116,32,98,117,105,108,116,32,119,105,116,104,32,108,105,98,101,120,112,97,116,46,32,84,97,98,108,101,32,102,111,114,109,97,116,116,105,110,103,32,105,115,32,110,111,116,32,97,118,97,105,108,97,98,108,101,46,10,0,0,0,0,47,114,101,100,115,53,47,52,0,0,0,0,0,0,0,0,47,114,101,100,115,53,47,51,0,0,0,0,0,0,0,0,47,98,114,98,103,56,47,53,0,0,0,0,0,0,0,0,47,114,101,100,115,53,47,50,0,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,60,68,73,86,32,105,100,61,39,95,86,77,76,49,95,39,32,115,116,121,108,101,61,34,112,111,115,105,116,105,111,110,58,114,101,108,97,116,105,118,101,59,32,100,105,115,112,108,97,121,58,105,110,108,105,110,101,59,32,118,105,115,105,98,105,108,105,116,121,58,104,105,100,100,101,110,0,0,0,0,115,101,97,103,114,101,101,110,0,0,0,0,0,0,0,0,9,9,9,40,92,40,41,32,115,104,111,119,32,105,32,115,116,114,32,99,118,115,32,115,104,111,119,32,40,44,41,32,115,104,111,119,32,106,32,115,116,114,32,99,118,115,32,115,104,111,119,32,40,92,41,41,32,115,104,111,119,0,0,0,47,47,42,42,42,32,112,111,108,121,108,105,110,101,10,0,47,114,101,100,115,53,47,49,0,0,0,0,0,0,0,0,47,114,101,100,115,52,47,52,0,0,0,0,0,0,0,0,99,114,97,114,114,0,0,0,47,114,101,100,115,52,47,51,0,0,0,0,0,0,0,0,73,83,79,56,56,53,57,45,49,0,0,0,0,0,0,0,83,97,116,117,114,100,97,121,0,0,0,0,0,0,0,0,47,114,101,100,115,52,47,50,0,0,0,0,0,0,0,0,47,114,101,100,115,52,47,49,0,0,0,0,0,0,0,0,47,114,101,100,115,51,47,51,0,0,0,0,0,0,0,0,37,99,37,108,100,0,0,0,47,114,101,100,115,51,47,50,0,0,0,0,0,0,0,0,47,114,101,100,115,51,47,49,0,0,0,0,0,0,0,0,47,114,100,121,108,103,110,57,47,57,0,0,0,0,0,0,70,114,105,100,97,121,0,0,47,114,100,121,108,103,110,57,47,56,0,0,0,0,0,0,60,66,79,68,89,32,111,110,108,111,97,100,61,39,98,114,111,119,115,101,114,99,104,101,99,107,40,41,59,39,62,10,0,0,0,0,0,0,0,0,115,97,110,100,121,98,114,111,119,110,0,0,0,0,0,0,9,9,9,48,32,48,32,109,111,118,101,116,111,0,0,0,47,47,42,42,42,32,99,111,109,109,101,110,116,58,32,37,115,10,0,0,0,0,0,0,47,114,100,121,108,103,110,57,47,55,0,0,0,0,0,0,47,97,99,99,101,110,116,53,47,52,0,0,0,0,0,0,47,114,100,121,108,103,110,57,47,54,0,0,0,0,0,0,99,111,112,121,0,0,0,0,47,114,100,121,108,103,110,57,47,53,0,0,0,0,0,0,73,83,79,95,56,56,53,57,45,49,0,0,0,0,0,0,47,114,100,121,108,103,110,57,47,52,0,0,0,0,0,0,47,114,100,121,108,103,110,57,47,51,0,0,0,0,0,0,47,114,100,121,108,103,110,57,47,50,0,0,0,0,0,0,47,98,114,98,103,56,47,52,0,0,0,0,0,0,0,0,100,111,117,98,108,101,0,0,47,114,100,121,108,103,110,57,47,49,0,0,0,0,0,0,47,114,100,121,108,103,110,56,47,56,0,0,0,0,0,0,98,108,117,101,0,0,0,0,47,114,100,121,108,103,110,56,47,55,0,0,0,0,0,0,84,104,117,114,115,100,97,121,0,0,0,0,0,0,0,0,47,114,100,121,108,103,110,56,47,54,0,0,0,0,0,0,60,47,72,69,65,68,62,0,115,97,108,109,111,110,0,0,9,9,9,99,111,111,114,100,102,111,110,116,32,115,101,116,102,111,110,116,0,0,0,0,47,114,100,121,108,103,110,56,47,53,0,0,0,0,0,0,98,108,117,101,118,105,111,108,101,116,0,0,0,0,0,0,118,105,101,119,66,111,120,0,47,114,100,121,108,103,110,56,47,52,0,0,0,0,0,0,99,111,110,103,0,0,0,0,47,114,100,121,108,103,110,56,47,51,0,0,0,0,0,0,108,49,0,0,0,0,0,0,47,114,100,121,108,103,110,56,47,50,0,0,0,0,0,0,47,114,100,121,108,103,110,56,47,49,0,0,0,0,0,0,98,114,111,110,122,101,50,0,47,114,100,121,108,103,110,55,47,55,0,0,0,0,0,0,47,98,114,98,103,56,47,51,0,0,0,0,0,0,0,0,67,111,117,108,100,32,110,111,116,32,112,97,114,115,101,32,34,95,100,114,97,119,95,34,32,97,116,116,114,105,98,117,116,101,32,105,110,32,103,114,97,112,104,32,37,115,10,0,46,112,115,32,37,100,42,92,110,40,83,70,117,47,37,46,48,102,117,10,0,0,0,0,46,46,47,46,46,47,46,46,47,112,108,117,103,105,110,47,99,111,114,101,47,103,118,114,101,110,100,101,114,95,99,111,114,101,95,109,97,112,46,99,0,0,0,0,0,0,0,0,47,114,100,121,108,103,110,55,47,54,0,0,0,0,0,0,47,114,100,121,108,103,110,55,47,53,0,0,0,0,0,0,106,112,101,103,58,118,114,109,108,0,0,0,0,0,0,0,47,114,100,121,108,103,110,55,47,52,0,0,0,0,0,0,87,101,100,110,101,115,100,97,121,0,0,0,0,0,0,0,47,114,100,121,108,103,110,55,47,51,0,0,0,0,0,0,32,32,32,60,47,83,67,82,73,80,84,62,10,0,0,0,115,97,100,100,108,101,98,114,111,119,110,0,0,0,0,0,9,9,103,115,97,118,101,0,47,114,100,121,108,103,110,55,47,50,0,0,0,0,0,0,47,114,100,121,108,103,110,55,47,49,0,0,0,0,0,0,98,32,61,61,32,110,0,0,99,108,117,98,115,0,0,0,115,101,108,102,0,0,0,0,47,114,100,121,108,103,110,54,47,54,0,0,0,0,0,0,108,97,116,105,110,49,0,0,47,114,100,121,108,103,110,54,47,53,0,0,0,0,0,0,47,114,100,121,108,103,110,54,47,52,0,0,0,0,0,0,75,80,95,83,117,98,116,114,97,99,116,0,0,0,0,0,105,32,61,61,32,100,101,103,0,0,0,0,0,0,0,0,47,114,100,121,108,103,110,54,47,51,0,0,0,0,0,0,47,98,114,98,103,56,47,50,0,0,0,0,0,0,0,0,114,97,110,107,105,110,103,58,32,102,97,105,108,117,114,101,32,116,111,32,99,114,101,97,116,101,32,115,116,114,111,110,103,32,99,111,110,115,116,114,97,105,110,116,32,101,100,103,101,32,98,101,116,119,101,101,110,32,110,111,100,101,115,32,37,115,32,97,110,100,32,37,115,10,0,0,0,0,0,0,111,114,100,101,114,105,110,103,32,39,37,115,39,32,110,111,116,32,114,101,99,111,103,110,105,122,101,100,32,102,111,114,32,110,111,100,101,32,39,37,115,39,46,10,0,0,0,0,47,114,100,121,108,103,110,54,47,50,0,0,0,0,0,0,110,111,110,97,109,101,46,103,118,0,0,0,0,0,0,0,37,108,102,44,37,100,0,0,40,37,46,53,103,44,37,46,53,103,41,0,0,0,0,0,47,114,100,121,108,103,110,54,47,49,0,0,0,0,0,0,47,114,100,121,108,103,110,53,47,53,0,0,0,0,0,0,99,97,110,110,111,116,32,114,101,97,108,108,111,99,32,112,110,108,112,115,0,0,0,0,84,117,101,115,100,97,121,0,47,114,100,121,108,103,110,53,47,52,0,0,0,0,0,0,32,32,32,125,10,0,0,0,114,111,121,97,108,98,108,117,101,0,0,0,0,0,0,0,9,110,112,97,103,101,115,32,49,32,103,116,32,123,0,0,47,114,100,121,108,103,110,53,47,51,0,0,0,0,0,0,37,115,32,105,115,32,98,105,103,103,101,114,32,116,104,97,110,32,73,32,99,97,110,32,104,97,110,100,108,101,46,10,0,0,0,0,0,0,0,0,47,114,100,121,108,103,110,53,47,50,0,0,0,0,0,0,99,105,114,99,0,0,0,0,47,114,100,121,108,103,110,53,47,49,0,0,0,0,0,0,108,97,116,105,110,45,49,0,47,114,100,121,108,103,110,52,47,52,0,0,0,0,0,0,37,37,69,110,100,68,111,99,117,109,101,110,116,10,0,0,47,114,100,121,108,103,110,52,47,51,0,0,0,0,0,0,114,101,99,116,115,0,0,0,47,114,100,121,108,103,110,52,47,50,0,0,0,0,0,0,47,98,114,98,103,56,47,49,0,0,0,0,0,0,0,0,47,114,100,121,108,103,110,52,47,49,0,0,0,0,0,0,47,114,100,121,108,103,110,51,47,51,0,0,0,0,0,0,47,114,100,121,108,103,110,51,47,50,0,0,0,0,0,0,77,111,110,100,97,121,0,0,47,114,100,121,108,103,110,51,47,49,0,0,0,0,0,0,32,32,32,32,32,125,10,0,114,111,115,121,98,114,111,119,110,0,0,0,0,0,0,0,9,47,115,116,114,32,49,48,32,115,116,114,105,110,103,32,100,101,102,0,0,0,0,0,47,114,100,121,108,103,110,49,49,47,57,0,0,0,0,0,47,114,100,121,108,103,110,49,49,47,56,0,0,0,0,0,102,97,116,97,108,32,102,108,101,120,32,115,99,97,110,110,101,114,32,105,110,116,101,114,110,97,108,32,101,114,114,111,114,45,45,101,110,100,32,111,102,32,98,117,102,102,101,114,32,109,105,115,115,101,100,0,99,104,105,0,0,0,0,0,47,114,100,121,108,103,110,49,49,47,55,0,0,0,0,0,117,116,102,45,56,0,0,0,47,114,100,121,108,103,110,49,49,47,54,0,0,0,0,0,47,114,100,121,108,103,110,49,49,47,53,0,0,0,0,0,47,114,100,121,108,103,110,49,49,47,52,0,0,0,0,0,47,98,114,98,103,55,47,55,0,0,0,0,0,0,0,0,47,114,100,121,108,103,110,49,49,47,51,0,0,0,0,0,47,114,100,121,108,103,110,49,49,47,50,0,0,0,0,0,47,114,100,121,108,103,110,49,49,47,49,49,0,0,0,0,83,117,110,100,97,121,0,0,47,114,100,121,108,103,110,49,49,47,49,48,0,0,0,0,32,32,32,32,32,125,101,108,115,101,123,10,0,0,0,0,9,47,105,32,101,120,99,104,32,100,101,102,0,0,0,0,118,105,111,108,101,116,114,101,100,0,0,0,0,0,0,0,47,114,100,121,108,103,110,49,49,47,49,0,0,0,0,0,47,114,100,121,108,103,110,49,48,47,57,0,0,0,0,0,99,101,110,116,0,0,0,0,47,114,100,121,108,103,110,49,48,47,56,0,0,0,0,0,99,104,97,114,115,101,116,0,112,97,103,101,37,100,44,37,100,95,0,0,0,0,0,0,47,114,100,121,108,103,110,49,48,47,55,0,0,0,0,0,47,114,100,121,108,103,110,49,48,47,54,0,0,0,0,0,47,114,100,121,108,103,110,49,48,47,53,0,0,0,0,0,47,98,114,98,103,55,47,54,0,0,0,0,0,0,0,0,47,114,100,121,108,103,110,49,48,47,52,0,0,0,0,0,47,114,100,121,108,103,110,49,48,47,51,0,0,0,0,0,47,114,100,121,108,103,110,49,48,47,50,0,0,0,0,0,58,32,0,0,0,0,0,0,122,0,0,0,0,0,0,0,47,114,100,121,108,103,110,49,48,47,49,48,0,0,0,0,32,32,32,32,32,32,32,32,32,32,32,105,116,101,109,46,115,116,121,108,101,46,118,105,115,105,98,105,108,105,116,121,61,39,104,105,100,100,101,110,39,59,10,0,0,0,0,0,9,47,106,32,101,120,99,104,32,100,101,102,0,0,0,0,47,114,100,121,108,103,110,49,48,47,49,0,0,0,0,0,47,114,100,121,108,98,117,57,47,57,0,0,0,0,0,0,99,101,100,105,108,0,0,0,47,114,100,121,108,98,117,57,47,56,0,0,0,0,0,0,102,105,108,108,0,0,0,0,73,108,108,101,103,97,108,32,108,101,110,103,116,104,32,118,97,108,117,101,32,105,110,32,34,37,115,34,32,99,111,108,111,114,32,97,116,116,114,105,98,117,116,101,32,0,0,0,47,114,100,121,108,98,117,57,47,55,0,0,0,0,0,0,47,114,100,121,108,98,117,57,47,54,0,0,0,0,0,0,47,114,100,121,108,98,117,57,47,53,0,0,0,0,0,0,47,98,114,98,103,55,47,53,0,0,0,0,0,0,0,0,103,118,82,101,110,100,101,114,83,116,114,105,110,103,0,0,47,114,100,121,108,98,117,57,47,52,0,0,0,0,0,0,47,114,100,121,108,98,117,57,47,51,0,0,0,0,0,0,47,114,100,121,108,98,117,57,47,50,0,0,0,0,0,0,97,116,116,114,105,98,117,116,101,32,109,97,99,114,111,115,32,110,111,116,32,105,109,112,108,101,109,101,110,116,101,100,0,0,0,0,0,0,0,0,83,0,0,0,97,0,0,0,116,0,0,0,0,0,0,0,47,114,100,121,108,98,117,57,47,49,0,0,0,0,0,0,32,32,32,32,32,32,32,32,32,105,116,101,109,32,61,32,100,111,99,117,109,101,110,116,46,103,101,116,69,108,101,109,101,110,116,66,121,73,100,40,86,77,76,110,111,91,120,93,41,59,10,0,0,0,0,0,47,114,100,121,108,98,117,56,47,56,0,0,0,0,0,0,112,111,119,100,101,114,98,108,117,101,0,0,0,0,0,0,9,47,110,112,97,103,101,115,32,101,120,99,104,32,100,101,102,0,0,0,0,0,0,0,118,101,114,121,95,108,105,103,104,116,95,112,117,114,112,108,101,0,0,0,0,0,0,0,112,101,110,100,32,100,105,99,116,111,102,32,97,32,98,97,100,32,111,98,106,101,99,116,0,0,0,0,0,0,0,0,47,114,100,121,108,98,117,56,47,55,0,0,0,0,0,0,99,99,101,100,105,108,0,0,47,114,100,121,108,98,117,56,47,54,0,0,0,0,0,0,101,120,112,97,110,100,0,0,84,111,116,97,108,32,115,105,122,101,32,62,32,49,32,105,110,32,34,37,115,34,32,99,111,108,111,114,32,115,112,101,99,32,0,0,0,0,0,0,47,114,100,121,108,98,117,56,47,53,0,0,0,0,0,0,47,114,100,121,108,98,117,56,47,52,0,0,0,0,0,0,47,114,100,121,108,98,117,56,47,51,0,0,0,0,0,0,47,98,114,98,103,55,47,52,0,0,0,0,0,0,0,0,47,114,100,121,108,98,117,56,47,50,0,0,0,0,0,0,47,114,100,121,108,98,117,56,47,49,0,0,0,0,0,0,47,114,100,121,108,98,117,55,47,55,0,0,0,0,0,0,47,114,100,121,108,98,117,55,47,54,0,0,0,0,0,0,32,32,32,32,32,32,32,102,111,114,32,40,120,32,105,110,32,86,77,76,110,111,41,123,10,0,0,0,0,0,0,0,112,108,117,109,0,0,0,0,47,98,101,103,105,110,112,97,103,101,32,123,9,37,32,105,32,106,32,110,112,97,103,101,115,0,0,0,0,0,0,0,118,101,114,121,100,97,114,107,98,114,111,119,110,0,0,0,47,114,100,121,108,98,117,55,47,53,0,0,0,0,0,0,47,114,100,121,108,98,117,55,47,52,0,0,0,0,0,0,99,97,112,0,0,0,0,0,47,114,100,121,108,98,117,55,47,51,0,0,0,0,0,0,70,0,0,0,114,0,0,0,105,0,0,0,0,0,0,0,47,114,100,121,108,98,117,55,47,50,0,0,0,0,0,0,47,114,100,121,108,98,117,55,47,49,0,0,0,0,0,0,47,114,100,121,108,98,117,54,47,54,0,0,0,0,0,0,47,98,114,98,103,55,47,51,0,0,0,0,0,0,0,0,47,114,100,121,108,98,117,54,47,53,0,0,0,0,0,0,97,103,100,101,108,101,116,101,32,111,110,32,98,97,100,32,111,98,106,101,99,116,0,0,47,114,100,121,108,98,117,54,47,52,0,0,0,0,0,0,47,114,100,121,108,98,117,54,47,51,0,0,0,0,0,0,84,0,0,0,104,0,0,0,117,0,0,0,0,0,0,0,47,114,100,121,108,98,117,54,47,50,0,0,0,0,0,0,32,32,32,32,32,32,32,125,10,0,0,0,0,0,0,0,112,105,110,107,0,0,0,0,47,110,111,112,99,111,108,111,114,32,123,112,111,112,32,112,111,112,32,112,111,112,125,32,98,105,110,100,32,100,101,102,0,0,0,0,0,0,0,0,47,114,100,121,108,98,117,54,47,49,0,0,0,0,0,0,47,114,100,121,108,98,117,53,47,53,0,0,0,0,0,0,98,117,108,108,0,0,0,0,47,114,100,121,108,98,117,53,47,52,0,0,0,0,0,0,97,117,116,111,0,0,0,0,108,97,121,101,114,115,32,110,111,116,32,115,117,112,112,111,114,116,101,100,32,105,110,32,37,115,32,111,117,116,112,117,116,10,0,0,0,0,0,0,47,114,100,121,108,98,117,53,47,51,0,0,0,0,0,0,47,114,100,121,108,98,117,53,47,50,0,0,0,0,0,0,47,114,100,121,108,98,117,53,47,49,0,0,0,0,0,0,47,98,114,98,103,55,47,50,0,0,0,0,0,0,0,0,47,114,100,121,108,98,117,52,47,52,0,0,0,0,0,0,47,114,100,121,108,98,117,52,47,51,0,0,0,0,0,0,47,114,100,121,108,98,117,52,47,50,0,0,0,0,0,0,87,0,0,0,101,0,0,0,100,0,0,0,0,0,0,0,47,114,100,121,108,98,117,52,47,49,0,0,0,0,0,0,32,32,32,32,32,32,32,32,32,125,10,0,0,0,0,0,112,101,114,117,0,0,0,0,47,103,114,97,112,104,99,111,108,111,114,32,123,32,115,101,116,104,115,98,99,111,108,111,114,32,125,32,98,105,110,100,32,100,101,102,0,0,0,0,47,114,100,121,108,98,117,51,47,51,0,0,0,0,0,0,47,97,99,99,101,110,116,53,47,51,0,0,0,0,0,0,47,114,100,121,108,98,117,51,47,50,0,0,0,0,0,0,98,114,118,98,97,114,0,0,47,114,100,121,108,98,117,51,47,49,0,0,0,0,0,0,114,97,116,105,111,0,0,0,73,109,97,103,101,115,32,117,110,115,117,112,112,111,114,116,101,100,32,105,110,32,34,98,97,99,107,103,114,111,117,110,100,34,32,97,116,116,114,105,98,117,116,101,10,0,0,0,99,117,114,118,101,0,0,0,47,114,100,121,108,98,117,49,49,47,57,0,0,0,0,0,47,114,100,121,108,98,117,49,49,47,56,0,0,0,0,0,47,114,100,121,108,98,117,49,49,47,55,0,0,0,0,0,47,98,114,98,103,55,47,49,0,0,0,0,0,0,0,0,102,108,111,97,116,0,0,0,47,114,100,121,108,98,117,49,49,47,54,0,0,0,0,0,47,114,100,121,108,98,117,49,49,47,53,0,0,0,0,0,98,108,97,99,107,0,0,0,101,100,103,101,0,0,0,0,47,114,100,121,108,98,117,49,49,47,52,0,0,0,0,0,84,0,0,0,117,0,0,0,101,0,0,0,0,0,0,0,47,114,100,121,108,98,117,49,49,47,51,0,0,0,0,0,32,32,32,32,32,32,32,32,32,32,32,105,116,101,109,46,115,116,121,108,101,46,118,105,115,105,98,105,108,105,116,121,61,39,118,105,115,105,98,108,101,39,59,10,0,0,0,0,112,101,97,99,104,112,117,102,102,0,0,0,0,0,0,0,47,101,100,103,101,99,111,108,111,114,32,123,32,115,101,116,104,115,98,99,111,108,111,114,32,125,32,98,105,110,100,32,100,101,102,0,0,0,0,0,47,114,100,121,108,98,117,49,49,47,50,0,0,0,0,0,47,114,100,121,108,98,117,49,49,47,49,49,0,0,0,0,47,114,100,121,108,98,117,49,49,47,49,48,0,0,0,0,98,101,116,97,0,0,0,0,37,108,102,37,99,0,0,0,117,110,115,112,101,99,105,102,105,101,100,32,105,111,115,116,114,101,97,109,95,99,97,116,101,103,111,114,121,32,101,114,114,111,114,0,0,0,0,0,105,110,118,105,115,0,0,0,47,114,100,121,108,98,117,49,49,47,49,0,0,0,0,0,103,118,114,101,110,100,101,114,95,115,101,116,95,115,116,121,108,101,58,32,117,110,115,117,112,112,111,114,116,101,100,32,115,116,121,108,101,32,37,115,32,45,32,105,103,110,111,114,105,110,103,10,0,0,0,0,47,114,100,121,108,98,117,49,48,47,57,0,0,0,0,0,98,114,111,110,122,101,0,0,47,114,100,121,108,98,117,49,48,47,56,0,0,0,0,0,47,98,114,98,103,54,47,54,0,0,0,0,0,0,0,0,46,102,116,32,37,115,10,0,47,114,100,121,108,98,117,49,48,47,55,0,0,0,0,0,37,100,32,37,100,32,37,100,32,37,100,32,37,100,32,37,100,32,37,100,32,37,100,32,37,100,32,37,46,49,102,32,37,100,32,37,100,32,37,100,32,37,100,32,37,100,32,37,100,10,0,0,0,0,0,0,32,37,100,32,0,0,0,0,47,114,100,121,108,98,117,49,48,47,54,0,0,0,0,0,103,105,102,58,118,114,109,108,0,0,0,0,0,0,0,0,47,114,100,121,108,98,117,49,48,47,53,0,0,0,0,0,77,0,0,0,111,0,0,0,110,0,0,0,0,0,0,0,47,114,100,121,108,98,117,49,48,47,52,0,0,0,0,0,32,32,32,32,32,32,32,32,32,105,102,32,40,105,116,101,109,41,32,123,10,0,0,0,112,97,112,97,121,97,119,104,105,112,0,0,0,0,0,0,47,110,111,100,101,99,111,108,111,114,32,123,32,115,101,116,104,115,98,99,111,108,111,114,32,125,32,98,105,110,100,32,100,101,102,0,0,0,0,0,115,117,109,109,101,114,115,107,121,0,0,0,0,0,0,0,47,114,100,121,108,98,117,49,48,47,51,0,0,0,0,0,47,114,100,121,108,98,117,49,48,47,50,0,0,0,0,0,98,100,113,117,111,0,0,0,115,116,97,114,116,0,0,0,47,114,100,121,108,98,117,49,48,47,49,48,0,0,0,0,47,114,100,121,108,98,117,49,48,47,49,0,0,0,0,0,47,114,100,112,117,57,47,57,0,0,0,0,0,0,0,0,105,100,120,32,61,61,32,115,122,0,0,0,0,0,0,0,47,114,100,112,117,57,47,56,0,0,0,0,0,0,0,0,47,98,114,98,103,54,47,53,0,0,0,0,0,0,0,0,95,119,101,97,107,95,37,100,0,0,0,0,0,0,0,0,111,114,100,101,114,105,110,103,32,39,37,115,39,32,110,111,116,32,114,101,99,111,103,110,105,122,101,100,46,10,0,0,47,114,100,112,117,57,47,55,0,0,0,0,0,0,0,0,46,37,100,0,0,0,0,0,97,117,120,103,0,0,0,0,97,115,112,101,99,116,0,0,47,114,100,112,117,57,47,54,0,0,0,0,0,0,0,0,47,114,100,112,117,57,47,53,0,0,0,0,0,0,0,0,99,97,110,110,111,116,32,114,101,97,108,108,111,99,32,112,110,108,115,0,0,0,0,0,83,0,0,0,117,0,0,0,110,0,0,0,0,0,0,0,47,114,100,112,117,57,47,52,0,0,0,0,0,0,0,0,32,32,32,32,32,32,32,32,32,105,116,101,109,32,61,32,100,111,99,117,109,101,110,116,46,103,101,116,69,108,101,109,101,110,116,66,121,73,100,40,86,77,76,121,101,115,91,120,93,41,59,10,0,0,0,0,112,97,108,101,118,105,111,108,101,116,114,101,100,0,0,0,37,32,104,111,111,107,115,32,102,111,114,32,115,101,116,116,105,110,103,32,99,111,108,111,114,32,0,0,0,0,0,0,47,114,100,112,117,57,47,51,0,0,0,0,0,0,0,0,47,114,100,112,117,57,47,50,0,0,0,0,0,0,0,0,97,117,109,108,0,0,0,0,101,108,108,105,112,115,101,0,47,114,100,112,117,57,47,49,0,0,0,0,0,0,0,0,108,97,98,101,108,106,117,115,116,0,0,0,0,0,0,0,47,114,100,112,117,56,47,56,0,0,0,0,0,0,0,0,37,37,66,101,103,105,110,68,111,99,117,109,101,110,116,58,10,0,0,0,0,0,0,0,47,114,100,112,117,56,47,55,0,0,0,0,0,0,0,0,47,114,100,112,117,56,47,54,0,0,0,0,0,0,0,0,47,98,114,98,103,54,47,52,0,0,0,0,0,0,0,0,47,114,100,112,117,56,47,53,0,0,0,0,0,0,0,0,95,37,108,100,95,83,85,83,80,69,67,84,0,0,0,0,47,114,100,112,117,56,47,52,0,0,0,0,0,0,0,0,47,114,100,112,117,56,47,51,0,0,0,0,0,0,0,0,115,97,110,115,45,83,101,114,105,102,0,0,0,0,0,0,83,0,0,0,97,0,0,0,116,0,0,0,117,0,0,0,114,0,0,0,100,0,0,0,97,0,0,0,121,0,0,0,0,0,0,0,0,0,0,0,47,114,100,112,117,56,47,50,0,0,0,0,0,0,0,0,32,32,32,32,32,32,32,102,111,114,32,40,120,32,105,110,32,86,77,76,121,101,115,41,123,10,0,0,0,0,0,0,112,97,108,101,116,117,114,113,117,111,105,115,101,0,0,0,47,116,97,112,101,114,101,100,32,123,32,125,32,98,105,110,100,32,100,101,102,0,0,0,47,114,100,112,117,56,47,49,0,0,0,0,0,0,0,0,47,114,100,112,117,55,47,55,0,0,0,0,0,0,0,0,97,116,105,108,100,101,0,0,111,117,116,32,111,102,32,100,121,110,97,109,105,99,32,109,101,109,111,114,121,32,105,110,32,97,97,103,101,110,115,117,114,101,95,98,117,102,102,101,114,95,115,116,97,99,107,40,41,0,0,0,0,0,0,0,47,114,100,112,117,55,47,54,0,0,0,0,0,0,0,0,47,114,100,112,117,55,47,53,0,0,0,0,0,0,0,0,47,114,100,112,117,55,47,52,0,0,0,0,0,0,0,0,47,114,100,112,117,55,47,51,0,0,0,0,0,0,0,0,47,98,114,98,103,54,47,51,0,0,0,0,0,0,0,0,47,114,100,112,117,55,47,50,0,0,0,0,0,0,0,0,47,114,100,112,117,55,47,49,0,0,0,0,0,0,0,0,47,114,100,112,117,54,47,54,0,0,0,0,0,0,0,0,70,0,0,0,114,0,0,0,105,0,0,0,100,0,0,0,97,0,0,0,121,0,0,0,0,0,0,0,0,0,0,0,47,114,100,112,117,54,47,53,0,0,0,0,0,0,0,0,32,32,32,32,32,32,105,102,32,40,105,101,118,101,114,115,62,61,53,41,123,10,0,0,112,97,108,101,103,114,101,101,110,0,0,0,0,0,0,0,47,100,105,97,103,111,110,97,108,115,32,123,32,125,32,98,105,110,100,32,100,101,102,0,115,112,105,99,121,112,105,110,107,0,0,0,0,0,0,0,47,114,100,112,117,54,47,52,0,0,0,0,0,0,0,0,47,114,100,112,117,54,47,51,0,0,0,0,0,0,0,0,116,105,108,100,101,0,0,0,97,115,121,109,112,0,0,0,47,114,100,112,117,54,47,50,0,0,0,0,0,0,0,0,47,114,100,112,117,54,47,49,0,0,0,0,0,0,0,0,47,114,100,112,117,53,47,53,0,0,0,0,0,0,0,0,47,114,100,112,117,53,47,52,0,0,0,0,0,0,0,0,47,98,114,98,103,54,47,50,0,0,0,0,0,0,0,0,47,114,100,112,117,53,47,51,0,0,0,0,0,0,0,0,47,114,100,112,117,53,47,50,0,0,0,0,0,0,0,0,47,114,100,112,117,53,47,49,0,0,0,0,0,0,0,0,84,0,0,0,104,0,0,0,117,0,0,0,114,0,0,0,115,0,0,0,100,0,0,0,97,0,0,0,121,0,0,0,0,0,0,0,0,0,0,0,47,114,100,112,117,52,47,52,0,0,0,0,0,0,0,0,32,32,32,32,32,32,125,10,0,0,0,0,0,0,0,0,112,97,108,101,103,111,108,100,101,110,114,111,100,0,0,0,47,114,111,117,110,100,101,100,32,123,32,125,32,98,105,110,100,32,100,101,102,0,0,0,47,114,100,112,117,52,47,51,0,0,0,0,0,0,0,0,47,114,100,112,117,52,47,50,0,0,0,0,0,0,0,0,97,114,105,110,103,0,0,0,47,114,100,112,117,52,47,49,0,0,0,0,0,0,0,0,47,114,100,112,117,51,47,51,0,0,0,0,0,0,0,0,47,114,100,112,117,51,47,50,0,0,0,0,0,0,0,0,47,114,100,112,117,51,47,49,0,0,0,0,0,0,0,0,47,98,114,98,103,54,47,49,0,0,0,0,0,0,0,0,103,118,76,97,121,111,117,116,0,0,0,0,0,0,0,0,47,114,100,103,121,57,47,57,0,0,0,0,0,0,0,0,47,114,100,103,121,57,47,56,0,0,0,0,0,0,0,0,47,114,100,103,121,57,47,55,0,0,0,0,0,0,0,0,47,114,100,103,121,57,47,54,0,0,0,0,0,0,0,0,87,0,0,0,101,0,0,0,100,0,0,0,110,0,0,0,101,0,0,0,115,0,0,0,100,0,0,0,97,0,0,0,121,0,0,0,0,0,0,0,37,115,58,37,115,0,0,0,32,32,32,32,32,32,32,32,32,105,101,118,101,114,115,61,32,112,97,114,115,101,73,110,116,32,40,117,97,46,115,117,98,115,116,114,105,110,103,32,40,109,115,105,101,43,53,44,32,117,97,46,105,110,100,101,120,79,102,32,40,39,46,39,44,32,109,115,105,101,32,41,41,41,10,0,0,0,0,0,111,114,99,104,105,100,0,0,47,117,110,102,105,108,108,101,100,32,123,32,125,32,98,105,110,100,32,100,101,102,0,0,47,114,100,103,121,57,47,53,0,0,0,0,0,0,0,0,47,114,100,103,121,57,47,52,0,0,0,0,0,0,0,0,116,97,98,0,0,0,0,0,97,110,103,0,0,0,0,0,47,114,100,103,121,57,47,51,0,0,0,0,0,0,0,0,47,114,100,103,121,57,47,50,0,0,0,0,0,0,0,0,47,114,100,103,121,57,47,49,0,0,0,0,0,0,0,0,47,114,100,103,121,56,47,56,0,0,0,0,0,0,0,0,47,98,114,98,103,53,47,53,0,0,0,0,0,0,0,0,47,114,100,103,121,56,47,55,0,0,0,0,0,0,0,0,47,114,100,103,121,56,47,54,0,0,0,0,0,0,0,0,47,114,100,103,121,56,47,53,0,0,0,0,0,0,0,0,84,0,0,0,117,0,0,0,101,0,0,0,115,0,0,0,100,0,0,0,97,0,0,0,121,0,0,0,0,0,0,0,47,114,100,103,121,56,47,52,0,0,0,0,0,0,0,0,32,32,32,32,32,32,105,102,32,40,32,109,115,105,101,32,62,32,48,32,41,123,32,32,32,32,32,32,47,47,32,73,102,32,73,110,116,101,114,110,101,116,32,69,120,112,108,111,114,101,114,44,32,114,101,116,117,114,110,32,118,101,114,115,105,111,110,32,110,117,109,98,101,114,10,0,0,0,0,0,111,114,97,110,103,101,114,101,100,0,0,0,0,0,0,0,47,102,105,108,108,101,100,32,123,32,125,32,98,105,110,100,32,100,101,102,0,0,0,0,47,114,100,103,121,56,47,51,0,0,0,0,0,0,0,0,47,114,100,103,121,56,47,50,0,0,0,0,0,0,0,0,97,110,100,0,0,0,0,0,47,114,100,103,121,56,47,49,0,0,0,0,0,0,0,0,47,114,100,103,121,55,47,55,0,0,0,0,0,0,0,0,47,114,100,103,121,55,47,54,0,0,0,0,0,0,0,0,47,114,100,103,121,55,47,53,0,0,0,0,0,0,0,0,47,98,114,98,103,53,47,52,0,0,0,0,0,0,0,0,47,114,100,103,121,55,47,52,0,0,0,0,0,0,0,0,47,114,100,103,121,55,47,51,0,0,0,0,0,0,0,0,47,114,100,103,121,55,47,50,0,0,0,0,0,0,0,0,77,0,0,0,111,0,0,0,110,0,0,0,100,0,0,0,97,0,0,0,121,0,0,0,0,0,0,0,0,0,0,0,47,114,100,103,121,55,47,49,0,0,0,0,0,0,0,0,32,32,32,32,32,32,118,97,114,32,86,77,76,110,111,61,110,101,119,32,65,114,114,97,121,40,39,95,110,111,116,86,77,76,49,95,39,44,39,95,110,111,116,86,77,76,50,95,39,41,59,10,0,0,0,0,111,114,97,110,103,101,0,0,47,98,111,108,100,32,123,32,50,32,115,101,116,108,105,110,101,119,105,100,116,104,32,125,32,98,105,110,100,32,100,101,102,0,0,0,0,0,0,0,47,114,100,103,121,54,47,54,0,0,0,0,0,0,0,0,47,114,100,103,121,54,47,53,0,0,0,0,0,0,0,0,97,109,112,0,0,0,0,0,47,114,100,103,121,54,47,52,0,0,0,0,0,0,0,0,73,83,79,45,56,56,53,57,45,49,0,0,0,0,0,0,104,101,97,100,116,111,111,108,116,105,112,0,0,0,0,0,47,114,100,103,121,54,47,51,0,0,0,0,0,0,0,0,47,114,100,103,121,54,47,50,0,0,0,0,0,0,0,0,47,114,100,103,121,54,47,49,0,0,0,0,0,0,0,0,47,98,114,98,103,53,47,51,0,0,0,0,0,0,0,0,47,114,100,103,121,53,47,53,0,0,0,0,0,0,0,0,47,114,100,103,121,53,47,52,0,0,0,0,0,0,0,0,47,114,100,103,121,53,47,51,0,0,0,0,0,0,0,0,47,98,114,98,103,53,47,50,0,0,0,0,0,0,0,0,83,0,0,0,117,0,0,0,110,0,0,0,100,0,0,0,97,0,0,0,121,0,0,0,0,0,0,0,0,0,0,0,47,114,100,103,121,53,47,50,0,0,0,0,0,0,0,0,32,32,32,32,32,32,118,97,114,32,86,77,76,121,101,115,61,110,101,119,32,65,114,114,97,121,40,39,95,86,77,76,49,95,39,44,39,95,86,77,76,50,95,39,41,59,10,0,111,108,105,118,101,100,114,97,98,0,0,0,0,0,0,0,47,105,110,118,105,115,32,123,47,102,105,108,108,32,123,110,101,119,112,97,116,104,125,32,100,101,102,32,47,115,116,114,111,107,101,32,123,110,101,119,112,97,116,104,125,32,100,101,102,32,47,115,104,111,119,32,123,112,111,112,32,110,101,119,112,97,116,104,125,32,100,101,102,125,32,98,105,110,100,32,100,101,102,0,0,0,0,0,115,101,109,105,83,119,101,101,116,67,104,111,99,0,0,0,47,114,100,103,121,53,47,49,0,0,0,0,0,0,0,0,47,97,99,99,101,110,116,53,47,50,0,0,0,0,0,0,47,114,100,103,121,52,47,52,0,0,0,0,0,0,0,0,47,114,100,103,121,52,47,51,0,0,0,0,0,0,0,0,116,97,105,108,116,111,111,108,116,105,112,0,0,0,0,0,109,112,116,121,0,0,0,0,47,114,100,103,121,52,47,50,0,0,0,0,0,0,0,0,47,114,100,103,121,52,47,49,0,0,0,0,0,0,0,0,47,114,100,103,121,51,47,51,0,0,0,0,0,0,0,0,117,110,115,105,103,110,101,100,32,108,111,110,103,0,0,0,47,114,100,103,121,51,47,50,0,0,0,0,0,0,0,0,47,114,100,103,121,51,47,49,0,0,0,0,0,0,0,0,97,113,117,97,0,0,0,0,47,114,100,103,121,49,49,47,57,0,0,0,0,0,0,0,47,114,100,103,121,49,49,47,56,0,0,0,0,0,0,0,32,32,32,32,32,32,118,97,114,32,105,116,101,109,59,10,0,0,0,0,0,0,0,0,47,100,111,116,116,101,100,32,123,32,91,49,32,73,110,118,83,99,97,108,101,70,97,99,116,111,114,32,109,117,108,32,54,32,73,110,118,83,99,97,108,101,70,97,99,116,111,114,32,109,117,108,93,32,48,32,115,101,116,100,97,115,104,32,125,32,98,105,110,100,32,100,101,102,0,0,0,0,0,0,47,114,100,103,121,49,49,47,55,0,0,0,0,0,0,0,98,108,97,110,99,104,101,100,97,108,109,111,110,100,0,0,112,116,0,0,0,0,0,0,47,114,100,103,121,49,49,47,54,0,0,0,0,0,0,0,97,108,101,102,115,121,109,0,47,114,100,103,121,49,49,47,53,0,0,0,0,0,0,0,108,97,98,101,108,116,111,111,108,116,105,112,0,0,0,0,37,115,32,0,0,0,0,0,47,114,100,103,121,49,49,47,52,0,0,0,0,0,0,0,47,114,100,103,121,49,49,47,51,0,0,0,0,0,0,0,98,114,105,103,104,116,103,111,108,100,0,0,0,0,0,0,47,114,100,103,121,49,49,47,50,0,0,0,0,0,0,0,47,98,114,98,103,53,47,49,0,0,0,0,0,0,0,0,32,37,100,44,37,100,0,0,47,114,100,103,121,49,49,47,49,49,0,0,0,0,0,0,35,32,37,115,10,0,0,0,120,100,111,116,58,120,100,111,116,0,0,0,0,0,0,0,47,114,100,103,121,49,49,47,49,48,0,0,0,0,0,0,112,110,103,58,118,114,109,108,0,0,0,0,0,0,0,0,47,114,100,103,121,49,49,47,49,0,0,0,0,0,0,0,108,97,98,101,108,95,115,99,104,101,109,101,32,61,32,37,100,32,62,32,52,32,58,32,105,103,110,111,114,105,110,103,10,0,0,0,0,0,0,0,68,101,99,0,0,0,0,0,40,33,106,99,110,41,32,38,38,32,40,33,118,97,108,41,0,0,0,0,0,0,0,0,47,114,100,103,121,49,48,47,57,0,0,0,0,0,0,0,32,32,32,32,32,32,118,97,114,32,105,101,118,101,114,115,59,10,0,0,0,0,0,0,111,108,100,108,97,99,101,0,47,100,97,115,104,101,100,32,123,32,91,57,32,73,110,118,83,99,97,108,101,70,97,99,116,111,114,32,109,117,108,32,100,117,112,32,93,32,48,32,115,101,116,100,97,115,104,32,125,32,98,105,110,100,32,100,101,102,0,0,0,0,0,0,115,99,97,114,108,101,116,0,47,114,100,103,121,49,48,47,56,0,0,0,0,0,0,0,47,114,100,103,121,49,48,47,55,0,0,0,0,0,0,0,97,103,114,97,118,101,0,0,110,111,100,101,32,39,37,115,39,44,32,103,114,97,112,104,32,39,37,115,39,32,115,105,122,101,32,116,111,111,32,115,109,97,108,108,32,102,111,114,32,108,97,98,101,108,10,0,120,108,112,0,0,0,0,0,47,114,100,103,121,49,48,47,54,0,0,0,0,0,0,0,105,100,0,0,0,0,0,0,101,100,103,101,116,111,111,108,116,105,112,0,0,0,0,0,47,114,100,103,121,49,48,47,53,0,0,0,0,0,0,0,47,114,100,103,121,49,48,47,52,0,0,0,0,0,0,0,101,115,101,112,0,0,0,0,46,46,47,46,46,47,46,46,47,108,105,98,47,102,100,112,103,101,110,47,108,97,121,111,117,116,46,99,0,0,0,0,75,80,95,65,100,100,0,0,47,114,100,103,121,49,48,47,51,0,0,0,0,0,0,0,47,98,114,98,103,52,47,52,0,0,0,0,0,0,0,0,99,111,109,112,97,99,116,0,105,110,0,0,0,0,0,0,47,114,100,103,121,49,48,47,50,0,0,0,0,0,0,0,109,101,109,111,114,121,32,97,108,108,111,99,97,116,105,111,110,32,102,97,105,108,117,114,101,10,0,0,0,0,0,0,108,97,98,101,108,0,0,0,115,101,103,109,101,110,116,32,91,37,115,44,37,115,93,32,100,111,101,115,32,110,111,116,32,105,110,116,101,114,115,101])
.concat([99,116,32,98,111,120,32,108,108,61,37,115,44,117,114,61,37,115,10,0,0,0,0,0,47,114,100,103,121,49,48,47,49,48,0,0,0,0,0,0,47,114,100,103,121,49,48,47,49,0,0,0,0,0,0,0,99,97,110,110,111,116,32,109,97,108,108,111,99,32,112,110,108,112,115,0,0,0,0,0,78,111,118,0,0,0,0,0,47,114,100,98,117,57,47,57,0,0,0,0,0,0,0,0,32,32,32,32,32,32,118,97,114,32,109,115,105,101,32,61,32,117,97,46,105,110,100,101,120,79,102,32,40,32,39,77,83,73,69,32,39,32,41,10,0,0,0,0,0,0,0,0,47,115,111,108,105,100,32,123,32,91,93,32,48,32,115,101,116,100,97,115,104,32,125,32,98,105,110,100,32,100,101,102,0,0,0,0,0,0,0,0,47,114,100,98,117,57,47,56,0,0,0,0,0,0,0,0,47,114,100,98,117,57,47,55,0,0,0,0,0,0,0,0,97,101,108,105,103,0,0,0,108,97,98,101,108,108,111,99,0,0,0,0,0,0,0,0,47,114,100,98,117,57,47,54,0,0,0,0,0,0,0,0,104,101,97,100,99,108,105,112,0,0,0,0,0,0,0,0,116,111,111,108,116,105,112,0,47,114,100,98,117,57,47,53,0,0,0,0,0,0,0,0,117,115,105,110,103,32,37,115,32,102,111,114,32,117,110,107,110,111,119,110,32,115,104,97,112,101,32,37,115,10,0,0,105,110,32,99,104,101,99,107,112,97,116,104,44,32,98,111,120,32,37,100,32,104,97,115,32,76,76,32,99,111,111,114,100,32,62,32,85,82,32,99,111,111,114,100,10,0,0,0,47,117,115,101,114,95,115,104,97,112,101,95,37,100,32,123,10,0,0,0,0,0,0,0,47,114,100,98,117,57,47,52,0,0,0,0,0,0,0,0,115,116,111,112,10,0,0,0,9,37,115,32,37,100,10,0,47,114,100,98,117,57,47,51,0,0,0,0,0,0,0,0,47,98,114,98,103,52,47,51,0,0,0,0,0,0,0,0,47,114,100,98,117,57,47,50,0,0,0,0,0,0,0,0,47,114,100,98,117,57,47,49,0,0,0,0,0,0,0,0,47,114,100,98,117,56,47,56,0,0,0,0,0,0,0,0,98,111,111,107,0,0,0,0,79,99,116,0,0,0,0,0,47,114,100,98,117,56,47,55,0,0,0,0,0,0,0,0,32,32,32,32,32,32,118,97,114,32,117,97,32,61,32,119,105,110,100,111,119,46,110,97,118,105,103,97,116,111,114,46,117,115,101,114,65,103,101,110,116,10,0,0,0,0,0,0,110,97,118,97,106,111,119,104,105,116,101,0,0,0,0,0,37,32,115,116,121,108,101,115,0,0,0,0,0,0,0,0,114,105,99,104,98,108,117,101,0,0,0,0,0,0,0,0,47,114,100,98,117,56,47,54,0,0,0,0,0,0,0,0,47,114,100,98,117,56,47,53,0,0,0,0,0,0,0,0,97,99,117,116,101,0,0,0,78,111,32,111,114,32,105,109,112,114,111,112,101,114,32,105,109,97,103,101,61,34,37,115,34,32,102,111,114,32,110,111,100,101,32,34,37,115,34,10,0,0,0,0,0,0,0,0,47,114,100,98,117,56,47,52,0,0,0,0,0,0,0,0,116,97,105,108,99,108,105,112,0,0,0,0,0,0,0,0,104,101,97,100,116,97,114,103,101,116,0,0,0,0,0,0,47,114,100,98,117,56,47,51,0,0,0,0,0,0,0,0,47,114,100,98,117,56,47,50,0,0,0,0,0,0,0,0,47,114,100,98,117,56,47,49,0,0,0,0,0,0,0,0,47,98,114,98,103,52,47,50,0,0,0,0,0,0,0,0,47,114,100,98,117,55,47,55,0,0,0,0,0,0,0,0,47,114,100,98,117,55,47,54,0,0,0,0,0,0,0,0,47,114,100,98,117,55,47,53,0,0,0,0,0,0,0,0,83,101,112,0,0,0,0,0,47,114,100,98,117,55,47,52,0,0,0,0,0,0,0,0,32,32,32,123,10,0,0,0,109,111,99,99,97,115,105,110,0,0,0,0,0,0,0,0,32,32,32,32,32,32,32,115,99,97,108,101,0,0,0,0,47,114,100,98,117,55,47,51,0,0,0,0,0,0,0,0,47,114,100,98,117,55,47,50,0,0,0,0,0,0,0,0,97,99,105,114,99,0,0,0,60,110,105,108,62,0,0,0,47,114,100,98,117,55,47,49,0,0,0,0,0,0,0,0,99,111,110,115,116,114,97,105,110,116,0,0,0,0,0,0,116,97,105,108,116,97,114,103,101,116,0,0,0,0,0,0,47,114,100,98,117,54,47,54,0,0,0,0,0,0,0,0,47,114,100,98,117,54,47,53,0,0,0,0,0,0,0,0,47,114,100,98,117,54,47,52,0,0,0,0,0,0,0,0,47,98,114,98,103,52,47,49,0,0,0,0,0,0,0,0,47,114,100,98,117,54,47,51,0,0,0,0,0,0,0,0,47,114,100,98,117,54,47,50,0,0,0,0,0,0,0,0,47,114,100,98,117,54,47,49,0,0,0,0,0,0,0,0,65,117,103,0,0,0,0,0,47,114,100,98,117,53,47,53,0,0,0,0,0,0,0,0,32,32,32,102,117,110,99,116,105,111,110,32,98,114,111,119,115,101,114,99,104,101,99,107,40,41,10,0,0,0,0,0,109,105,115,116,121,114,111,115,101,0,0,0,0,0,0,0,32,32,32,32,32,32,32,100,117,112,32,49,32,101,120,99,104,32,100,105,118,32,47,73,110,118,83,99,97,108,101,70,97,99,116,111,114,32,101,120,99,104,32,100,101,102,0,0,113,117,97,114,116,122,0,0,47,114,100,98,117,53,47,52,0,0,0,0,0,0,0,0,47,114,100,98,117,53,47,51,0,0,0,0,0,0,0,0,97,97,99,117,116,101,0,0,78,111,32,111,114,32,105,109,112,114,111,112,101,114,32,115,104,97,112,101,102,105,108,101,61,34,37,115,34,32,102,111,114,32,110,111,100,101,32,34,37,115,34,10,0,0,0,0,47,114,100,98,117,53,47,50,0,0,0,0,0,0,0,0,97,114,114,111,119,115,105,122,101,0,0,0,0,0,0,0,108,97,98,101,108,116,97,114,103,101,116,0,0,0,0,0,47,114,100,98,117,53,47,49,0,0,0,0,0,0,0,0,47,114,100,98,117,52,47,52,0,0,0,0,0,0,0,0,47,114,100,98,117,52,47,51,0,0,0,0,0,0,0,0,47,98,114,98,103,51,47,51,0,0,0,0,0,0,0,0,103,118,67,111,110,116,101,120,116,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,216,131,1,0,160,149,2,0,72,1,2,0,112,176,1,0,216,250,1,0,152,254,1,0,0,0,0,0,0,0,0,0,118,109,108,95,116,101,120,116,112,97,114,97,0,0,0,0,118,109,108,95,112,114,105,110,116,95,99,111,108,111,114,0,117,110,105,102,111,114,109,95,115,116,114,101,115,115,0,0,116,114,97,110,115,112,111,115,101,95,115,116,101,112,0,0,116,107,103,101,110,95,112,114,105,110,116,95,116,97,103,115,0,0,0,0,0,0,0,0,116,107,103,101,110,95,112,114,105,110,116,95,99,111,108,111,114,0,0,0,0,0,0,0,115,118,103,95,116,101,120,116,112,97,114,97,0,0,0,0,115,118,103,95,112,114,105,110,116,95,99,111,108,111,114,0,115,101,116,98,111,117,110,100,115,0,0,0,0,0,0,0,115,101,116,69,100,103,101,76,97,98,101,108,80,111,115,0,114,111,117,110,100,95,99,111,114,110,101,114,115,0,0,0,114,101,109,111,118,101,95,102,114,111,109,95,114,97,110,107,0,0,0,0,0,0,0,0,114,101,109,111,118,101,68,101,103,108,105,115,116,0,0,0,112,111,118,95,99,111,108,111,114,95,97,115,95,115,116,114,0,0,0,0,0,0,0,0,112,111,115,116,111,114,100,101,114,0,0,0,0,0,0,0,112,111,112,95,111,98,106,95,115,116,97,116,101,0,0,0,112,111,112,0,0,0,0,0,112,111,108,121,108,105,110,101,77,105,100,112,111,105,110,116,0,0,0,0,0,0,0,0,111,118,101,114,108,97,112,95,98,101,122,105,101,114,0,0,110,101,105,103,104,98,111,114,0,0,0,0,0,0,0,0,110,101,97,116,111,95,101,110,113,117,101,117,101,0,0,0,109,107,78,67,111,110,115,116,114,97,105,110,116,71,0,0,109,105,110,109,97,120,95,101,100,103,101,115,0,0,0,0,109,101,114,103,101,118,105,114,116,117,97,108,0,0,0,0,109,101,114,103,101,95,111,110,101,119,97,121,0,0,0,0,109,101,114,103,101,95,99,104,97,105,110,0,0,0,0,0,109,97,112,95,112,97,116,104,0,0,0,0,0,0,0,0,109,97,112,95,111,117,116,112,117,116,95,115,104,97,112,101,0,0,0,0,0,0,0,0,109,97,112,78,0,0,0,0,109,97,107,101,95,108,97,98,101,108,0,0,0,0,0,0,109,97,107,101,95,99,104,97,105,110,0,0,0,0,0,0,109,97,107,101,95,98,97,114,114,105,101,114,115,0,0,0,109,97,107,101,83,101,108,102,69,100,103,101,0,0,0,0,109,97,107,101,71,114,97,112,104,68,97,116,97,0,0,0,109,97,107,101,67,111,109,112,111,117,110,100,69,100,103,101,0,0,0,0,0,0,0,0,105,110,115,116,97,108,108,95,105,110,95,114,97,110,107,0,105,110,115,101,114,116,78,111,100,101,108,105,115,116,0,0,105,110,105,116,95,115,112,108,105,110,101,115,95,98,98,0,105,100,101,97,108,95,100,105,115,116,97,110,99,101,95,109,97,116,114,105,120,0,0,0,103,118,117,115,101,114,115,104,97,112,101,95,102,105,108,101,95,97,99,99,101,115,115,0,103,101,116,95,101,100,103,101,95,108,97,98,101,108,95,109,97,116,114,105,120,0,0,0,103,101,116,69,100,103,101,76,105,115,116,0,0,0,0,0,102,108,97,116,95,115,101,97,114,99,104,0,0,0,0,0,102,108,97,116,95,114,101,111,114,100,101,114,0,0,0,0,102,105,110,100,67,67,111,109,112,0,0,0,0,0,0,0,102,105,103,95,114,101,115,111,108,118,101,95,99,111,108,111,114,0,0,0,0,0,0,0,102,105,103,95,98,101,122,105,101,114,0,0,0,0,0,0,102,97,115,116,95,110,111,100,101,97,112,112,0,0,0,0,102,97,115,116,95,110,111,100,101,0,0,0,0,0,0,0,101,120,112,97,110,100,67,108,117,115,116,101,114,0,0,0,101,110,100,112,97,116,104,0,101,109,105,116,95,101,100,103,101,95,108,97,98,101,108,0,100,111,116,95,112,111,115,105,116,105,111,110,0,0,0,0,100,102,115,67,121,99,108,101,0,0,0,0,0,0,0,0,100,101,108,101,116,101,95,102,108,97,116,95,101,100,103,101,0,0,0,0,0,0,0,0,100,101,108,101,116,101,95,102,97,115,116,95,110,111,100,101,0,0,0,0,0,0,0,0,100,101,108,101,116,101,95,102,97,115,116,95,101,100,103,101,0,0,0,0,0,0,0,0,99,111,114,101,95,108,111,97,100,105,109,97,103,101,95,118,114,109,108,0,0,0,0,0,99,111,114,101,95,108,111,97,100,105,109,97,103,101,95,115,118,103,0,0,0,0,0,0,99,111,114,101,95,108,111,97,100,105,109,97,103,101,95,112,115,108,105,98,0,0,0,0,99,111,114,101,95,108,111,97,100,105,109,97,103,101,95,112,115,0,0,0,0,0,0,0,99,111,114,101,95,108,111,97,100,105,109,97,103,101,95,102,105,103,0,0,0,0,0,0,99,111,110,110,101,99,116,71,114,97,112,104,0,0,0,0,99,111,109,112,117,116,101,83,99,97,108,101,88,89,0,0,99,108,117,115,116,101,114,95,108,101,97,100,101,114,0,0,98,111,120,73,110,116,101,114,115,101,99,116,102,0,0,0,98,101,122,105,101,114,95,98,98,0,0,0,0,0,0,0,98,101,103,105,110,112,97,116,104,0,0,0,0,0,0,0,98,97,108,97,110,99,101,0,97,98,111,109,105,110,97,116,105,111,110,0,0,0,0,0,95,110,101,97,116,111,95,115,101,116,95,97,115,112,101,99,116,0,0,0,0,0,0,0,95,100,111,116,95,115,112,108,105,110,101,115,0,0,0,0,85,110,105,102,111,114,109,83,116,114,101,115,115,83,109,111,111,116,104,101,114,95,110,101,119,0,0,0,0,0,0,0,85,70,95,115,101,116,110,97,109,101,0,0,0,0,0,0,84,114,105,97,110,103,108,101,83,109,111,111,116,104,101,114,95,110,101,119,0,0,0,0,83,116,114,101,115,115,77,97,106,111,114,105,122,97,116,105,111,110,83,109,111,111,116,104,101,114,95,115,109,111,111,116,104,0,0,0,0,0,0,0,83,116,114,101,115,115,77,97,106,111,114,105,122,97,116,105,111,110,83,109,111,111,116,104,101,114,50,95,110,101,119,0,83,112,114,105,110,103,83,109,111,111,116,104,101,114,95,115,109,111,111,116,104,0,0,0,83,112,114,105,110,103,83,109,111,111,116,104,101,114,95,110,101,119,0,0,0,0,0,0,83,112,97,114,115,101,83,116,114,101,115,115,77,97,106,111,114,105,122,97,116,105,111,110,83,109,111,111,116,104,101,114,95,110,101,119,0,0,0,0,81,117,97,100,84,114,101,101,95,114,101,112,117,108,115,105,118,101,95,102,111,114,99,101,95,105,110,116,101,114,97,99,116,0,0,0,0,0,0,0,81,117,97,100,84,114,101,101,95,114,101,112,117,108,115,105,118,101,95,102,111,114,99,101,95,97,99,99,117,109,117,108,97,116,101,0,0,0,0,0,81,117,97,100,84,114,101,101,95,110,101,119,0,0,0,0,81,117,97,100,84,114,101,101,95,97,100,100,95,105,110,116,101,114,110,97,108,0,0,0,80,114,105,111,114,105,116,121,81,117,101,117,101,95,112,117,115,104,0,0,0,0,0,0,80,111,98,115,112,97,116,104,0,0,0,0,0,0,0,0,73,77,68,83,95,103,105,118,101,110,95,100,105,109,0,0,84,89,80,69,73,68,32,101,109,115,99,114,105,112,116,101,110,58,58,105,110,116,101,114,110,97,108,58,58,103,101,116,65,99,116,117,97,108,84,121,112,101,40,84,32,42,41,32,91,84,32,61,32,65,103,114,97,112,104,95,115,93,0,0,84,89,80,69,73,68,32,101,109,115,99,114,105,112,116,101,110,58,58,105,110,116,101,114,110,97,108,58,58,103,101,116,65,99,116,117,97,108,84,121,112,101,40,84,32,42,41,32,91,84,32,61,32,65,103,110,111,100,101,95,115,93,0,0,84,89,80,69,73,68,32,101,109,115,99,114,105,112,116,101,110,58,58,105,110,116,101,114,110,97,108,58,58,103,101,116,65,99,116,117,97,108,84,121,112,101,40,84,32,42,41,32,91,84,32,61,32,65,103,101,100,103,101,95,115,93,0,0,84,89,80,69,73,68,32,101,109,115,99,114,105,112,116,101,110,58,58,105,110,116,101,114,110,97,108,58,58,103,101,116,65,99,116,117,97,108,84,121,112,101,40,84,32,42,41,32,91,84,32,61,32,71,86,67,95,115,93,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,48,49,50,51,52,53,54,55,56,57,0,0,0,0,0,0,48,49,50,51,52,53,54,55,56,57,0,0,0,0,0,0,37,0,0,0,89,0,0,0,45,0,0,0,37,0,0,0,109,0,0,0,45,0,0,0,37,0,0,0,100,0,0,0,37,0,0,0,72,0,0,0,58,0,0,0,37,0,0,0,77,0,0,0,58,0,0,0,37,0,0,0,83,0,0,0,37,0,0,0,72,0,0,0,58,0,0,0,37,0,0,0,77,0,0,0,0,0,0,0,37,0,0,0,73,0,0,0,58,0,0,0,37,0,0,0,77,0,0,0,58,0,0,0,37,0,0,0,83,0,0,0,32,0,0,0,37,0,0,0,112,0,0,0,0,0,0,0,37,0,0,0,109,0,0,0,47,0,0,0,37,0,0,0,100,0,0,0,47,0,0,0,37,0,0,0,121,0,0,0,37,0,0,0,72,0,0,0,58,0,0,0,37,0,0,0,77,0,0,0,58,0,0,0,37,0,0,0,83,0,0,0,37,72,58,37,77,58,37,83,37,72,58,37,77,0,0,0,37,73,58,37,77,58,37,83,32,37,112,0,0,0,0,0,37,89,45,37,109,45,37,100,37,109,47,37,100,47,37,121,37,72,58,37,77,58,37,83,37,0,0,0,0,0,0,0,37,112,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,168,218,2,0,150,0,0,0,128,2,0,0,98,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,184,218,2,0,34,4,0,0,92,3,0,0,176,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,200,218,2,0,148,1,0,0,188,5,0,0,198,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,216,218,2,0,80,0,0,0,110,0,0,0,168,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,232,218,2,0,80,0,0,0,44,0,0,0,168,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,248,218,2,0,80,0,0,0,174,4,0,0,168,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,160,219,2,0,106,3,0,0,206,1,0,0,246,0,0,0,118,5,0,0,100,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,192,219,2,0,230,0,0,0,188,3,0,0,246,0,0,0,96,5,0,0,212,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,224,219,2,0,90,3,0,0,116,0,0,0,246,0,0,0,172,3,0,0,136,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,220,2,0,74,3,0,0,4,3,0,0,246,0,0,0,148,3,0,0,54,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,160,220,2,0,146,5,0,0,252,1,0,0,246,0,0,0,254,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,192,220,2,0,84,3,0,0,100,2,0,0,246,0,0,0,142,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,224,220,2,0,190,0,0,0,102,2,0,0,246,0,0,0,44,1,0,0,56,0,0,0,202,3,0,0,68,0,0,0,172,1,0,0,254,4,0,0,222,1,0,0,248,255,255,255,224,220,2,0,238,0,0,0,92,0,0,0,130,1,0,0,172,0,0,0,26,0,0,0,108,1,0,0,50,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,221,2,0,124,5,0,0,20,5,0,0,246,0,0,0,232,0,0,0,36,0,0,0,58,5,0,0,16,3,0,0,70,3,0,0,40,0,0,0,94,0,0,0,248,255,255,255,8,221,2,0,234,2,0,0,124,4,0,0,216,4,0,0,30,5,0,0,148,2,0,0,246,1,0,0,48,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,48,221,2,0,186,1,0,0,214,3,0,0,246,0,0,0,20,2,0,0,216,1,0,0,72,5,0,0,240,2,0,0,124,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,64,221,2,0,78,1,0,0,120,1,0,0,246,0,0,0,228,1,0,0,180,3,0,0,92,1,0,0,154,3,0,0,12,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,80,221,2,0,130,5,0,0,2,0,0,0,246,0,0,0,36,3,0,0,158,5,0,0,86,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,112,221,2,0,234,0,0,0,28,0,0,0,246,0,0,0,12,5,0,0,182,1,0,0,142,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,144,221,2,0,238,4,0,0,120,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,152,221,2,0,138,0,0,0,0,3,0,0,198,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,168,221,2,0,8,1,0,0,114,3,0,0,246,0,0,0,202,0,0,0,46,2,0,0,176,0,0,0,184,0,0,0,174,0,0,0,196,0,0,0,194,0,0,0,74,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,200,221,2,0,14,2,0,0,2,2,0,0,246,0,0,0,248,4,0,0,32,4,0,0,18,4,0,0,50,2,0,0,16,4,0,0,24,4,0,0,22,4,0,0,120,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,232,221,2,0,192,0,0,0,98,0,0,0,246,0,0,0,144,4,0,0,64,4,0,0,218,3,0,0,20,4,0,0,174,3,0,0,134,4,0,0,114,4,0,0,152,4,0,0,148,4,0,0,146,4,0,0,202,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,222,2,0,16,1,0,0,10,0,0,0,246,0,0,0,116,3,0,0,94,5,0,0,88,5,0,0,90,5,0,0,44,5,0,0,92,5,0,0,86,5,0,0,98,3,0,0,102,5,0,0,100,5,0,0,128,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,40,222,2,0,170,1,0,0,248,1,0,0,246,0,0,0,190,2,0,0,14,4,0,0,104,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,56,222,2,0,136,0,0,0,126,3,0,0,246,0,0,0,252,3,0,0,194,4,0,0,86,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,72,222,2,0,32,5,0,0,220,2,0,0,246,0,0,0,4,4,0,0,34,1,0,0,248,3,0,0,200,0,0,0,254,3,0,0,226,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,104,222,2,0,142,3,0,0,88,1,0,0,246,0,0,0,112,0,0,0,238,2,0,0,82,2,0,0,154,4,0,0,84,4,0,0,150,3,0,0,74,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,136,222,2,0,142,3,0,0,230,2,0,0,246,0,0,0,154,5,0,0,38,1,0,0,162,0,0,0,166,5,0,0,192,4,0,0,70,0,0,0,172,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,168,222,2,0,142,3,0,0,42,3,0,0,246,0,0,0,172,2,0,0,178,2,0,0,52,4,0,0,156,1,0,0,18,3,0,0,46,1,0,0,48,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,200,222,2,0,142,3,0,0,166,0,0,0,246,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,216,222,2,0,52,1,0,0,60,3,0,0,246,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,232,222,2,0,142,3,0,0,100,4,0,0,246,0,0,0,8,3,0,0,146,1,0,0,168,2,0,0,144,5,0,0,150,1,0,0,56,4,0,0,226,3,0,0,126,0,0,0,252,0,0,0,220,4,0,0,52,2,0,0,154,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,223,2,0,182,5,0,0,178,0,0,0,246,0,0,0,62,0,0,0,20,3,0,0,210,2,0,0,212,4,0,0,40,1,0,0,218,2,0,0,94,3,0,0,14,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,64,223,2,0,128,1,0,0,12,1,0,0,66,3,0,0,66,4,0,0,154,2,0,0,164,4,0,0,170,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,80,223,2,0,142,3,0,0,212,1,0,0,246,0,0,0,172,2,0,0,178,2,0,0,52,4,0,0,156,1,0,0,18,3,0,0,46,1,0,0,48,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,96,223,2,0,142,3,0,0,104,5,0,0,246,0,0,0,172,2,0,0,178,2,0,0,52,4,0,0,156,1,0,0,18,3,0,0,46,1,0,0,48,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,112,223,2,0,184,2,0,0,66,5,0,0,162,1,0,0,34,3,0,0,32,3,0,0,160,3,0,0,96,1,0,0,106,4,0,0,196,4,0,0,56,1,0,0,18,1,0,0,168,5,0,0,176,5,0,0,206,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,120,223,2,0,52,0,0,0,130,2,0,0,224,3,0,0,34,5,0,0,28,5,0,0,108,2,0,0,26,2,0,0,196,3,0,0,200,2,0,0,74,0,0,0,128,0,0,0,74,5,0,0,160,2,0,0,84,1,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,168,223,2,0,220,0,0,0,188,4,0,0,252,255,255,255,252,255,255,255,168,223,2,0,132,1,0,0,182,2,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,192,223,2,0,240,4,0,0,76,5,0,0,252,255,255,255,252,255,255,255,192,223,2,0,94,2,0,0,58,4,0,0,0,0,0,0,0,0,0,0,8,0,0,0,0,0,0,0,216,223,2,0,224,1,0,0,190,5,0,0,248,255,255,255,248,255,255,255,216,223,2,0,144,3,0,0,64,5,0,0,0,0,0,0,0,0,0,0,8,0,0,0,0,0,0,0,240,223,2,0,86,2,0,0,190,3,0,0,248,255,255,255,248,255,255,255,240,223,2,0,246,2,0,0,14,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,224,2,0,96,4,0,0,146,3,0,0,198,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,96,224,2,0,132,5,0,0,26,5,0,0,136,1,0,0,34,3,0,0,32,3,0,0,160,3,0,0,68,2,0,0,106,4,0,0,196,4,0,0,56,1,0,0,18,1,0,0,168,5,0,0,176,5,0,0,84,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,112,224,2,0,122,2,0,0,136,3,0,0,110,2,0,0,34,5,0,0,28,5,0,0,108,2,0,0,186,0,0,0,196,3,0,0,200,2,0,0,74,0,0,0,128,0,0,0,74,5,0,0,160,2,0,0,122,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,160,224,2,0,4,5,0,0,10,3,0,0,246,0,0,0,226,2,0,0,222,4,0,0,72,3,0,0,180,4,0,0,122,0,0,0,152,1,0,0,40,2,0,0,44,4,0,0,212,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,192,224,2,0,74,2,0,0,42,1,0,0,246,0,0,0,190,4,0,0,208,4,0,0,102,4,0,0,8,5,0,0,38,5,0,0,240,1,0,0,202,4,0,0,128,3,0,0,22,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,224,224,2,0,22,5,0,0,124,2,0,0,246,0,0,0,198,0,0,0,116,2,0,0,160,5,0,0,56,3,0,0,60,5,0,0,132,3,0,0,48,4,0,0,158,3,0,0,110,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,225,2,0,180,1,0,0,110,3,0,0,246,0,0,0,108,4,0,0,176,4,0,0,34,2,0,0,218,4,0,0,250,1,0,0,164,1,0,0,68,3,0,0,182,4,0,0,166,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,56,225,2,0,8,4,0,0,76,0,0,0,182,3,0,0,34,3,0,0,32,3,0,0,160,3,0,0,96,1,0,0,106,4,0,0,196,4,0,0,254,2,0,0,152,3,0,0,114,1,0,0,176,5,0,0,206,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,72,225,2,0,104,1,0,0,244,4,0,0,94,4,0,0,34,5,0,0,28,5,0,0,108,2,0,0,26,2,0,0,196,3,0,0,200,2,0,0,62,4,0,0,2,1,0,0,72,0,0,0,160,2,0,0,84,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,104,225,2,0,142,5,0,0,54,2,0,0,16,2,0,0,8,2,0,0,104,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,120,225,2,0,142,5,0,0,38,4,0,0,16,2,0,0,8,2,0,0,174,1,0,0,146,0,0,0,226,4,0,0,44,2,0,0,0,0,0,0,0,0,0,0,120,0,0,0,0,0,0,0,119,0,0,0,0,0,0,0,118,0,0,0,0,0,0,0,116,0,0,0,0,0,0,0,115,0,0,0,0,0,0,0,109,0,0,0,0,0,0,0,108,0,0,0,0,0,0,0,106,0,0,0,0,0,0,0,105,0,0,0,0,0,0,0,104,0,0,0,0,0,0,0,102,0,0,0,0,0,0,0,101,0,0,0,0,0,0,0,100,0,0,0,0,0,0,0,99,0,0,0,0,0,0,0,98,0,0,0,0,0,0,0,97,0,0,0,0,0,0,0,83,116,57,116,121,112,101,95,105,110,102,111,0,0,0,0,83,116,57,101,120,99,101,112,116,105,111,110,0,0,0,0,83,116,57,98,97,100,95,97,108,108,111,99,0,0,0,0,83,116,56,98,97,100,95,99,97,115,116,0,0,0,0,0,83,116,49,51,114,117,110,116,105,109,101,95,101,114,114,111,114,0,0,0,0,0,0,0,83,116,49,50,111,117,116,95,111,102,95,114,97,110,103,101,0,0,0,0,0,0,0,0,83,116,49,50,108,101,110,103,116,104,95,101,114,114,111,114,0,0,0,0,0,0,0,0,83,116,49,49,108,111,103,105,99,95,101,114,114,111,114,0,80,105,0,0,0,0,0,0,80,99,0,0,0,0,0,0,80,75,56,65,103,114,97,112,104,95,115,0,0,0,0,0,80,75,56,65,103,110,111,100,101,95,115,0,0,0,0,0,80,75,56,65,103,101,100,103,101,95,115,0,0,0,0,0,80,75,53,71,86,67,95,115,0,0,0,0,0,0,0,0,80,56,65,103,114,97,112,104,95,115,0,0,0,0,0,0,80,56,65,103,110,111,100,101,95,115,0,0,0,0,0,0,80,56,65,103,101,100,103,101,95,115,0,0,0,0,0,0,80,53,71,86,67,95,115,0,78,83,116,51,95,95,49,57,116,105,109,101,95,98,97,115,101,69,0,0,0,0,0,0,78,83,116,51,95,95,49,57,109,111,110,101,121,95,112,117,116,73,119,78,83,95,49,57,111,115,116,114,101,97,109,98,117,102,95,105,116,101,114,97,116,111,114,73,119,78,83,95,49,49,99,104,97,114,95,116,114,97,105,116,115,73,119,69,69,69,69,69,69,0,0,0,78,83,116,51,95,95,49,57,109,111,110,101,121,95,112,117,116,73,99,78,83,95,49,57,111,115,116,114,101,97,109,98,117,102,95,105,116,101,114,97,116,111,114,73,99,78,83,95,49,49,99,104,97,114,95,116,114,97,105,116,115,73,99,69,69,69,69,69,69,0,0,0,78,83,116,51,95,95,49,57,109,111,110,101,121,95,103,101,116,73,119,78,83,95,49,57,105,115,116,114,101,97,109,98,117,102,95,105,116,101,114,97,116,111,114,73,119,78,83,95,49,49,99,104,97,114,95,116,114,97,105,116,115,73,119,69,69,69,69,69,69,0,0,0,78,83,116,51,95,95,49,57,109,111,110,101,121,95,103,101,116,73,99,78,83,95,49,57,105,115,116,114,101,97,109,98,117,102,95,105,116,101,114,97,116,111,114,73,99,78,83,95,49,49,99,104,97,114,95,116,114,97,105,116,115,73,99,69,69,69,69,69,69,0,0,0,78,83,116,51,95,95,49,57,98,97,115,105,99,95,105,111,115,73,119,78,83,95,49,49,99,104,97,114,95,116,114,97,105,116,115,73,119,69,69,69,69,0,0,0,0,0,0,0,78,83,116,51,95,95,49,57,98,97,115,105,99,95,105,111,115,73,99,78,83,95,49,49,99,104,97,114,95,116,114,97,105,116,115,73,99,69,69,69,69,0,0,0,0,0,0,0,78,83,116,51,95,95,49,57,95,95,110,117,109,95,112,117,116,73,119,69,69,0,0,0,78,83,116,51,95,95,49,57,95,95,110,117,109,95,112,117,116,73,99,69,69,0,0,0,78,83,116,51,95,95,49,57,95,95,110,117,109,95,103,101,116,73,119,69,69,0,0,0,78,83,116,51,95,95,49,57,95,95,110,117,109,95,103,101,116,73,99,69,69,0,0,0,78,83,116,51,95,95,49,56,116,105,109,101,95,112,117,116,73,119,78,83,95,49,57,111,115,116,114,101,97,109,98,117,102,95,105,116,101,114,97,116,111,114,73,119,78,83,95,49,49,99,104,97,114,95,116,114,97,105,116,115,73,119,69,69,69,69,69,69,0,0,0,0,78,83,116,51,95,95,49,56,116,105,109,101,95,112,117,116])
.concat([73,99,78,83,95,49,57,111,115,116,114,101,97,109,98,117,102,95,105,116,101,114,97,116,111,114,73,99,78,83,95,49,49,99,104,97,114,95,116,114,97,105,116,115,73,99,69,69,69,69,69,69,0,0,0,0,78,83,116,51,95,95,49,56,116,105,109,101,95,103,101,116,73,119,78,83,95,49,57,105,115,116,114,101,97,109,98,117,102,95,105,116,101,114,97,116,111,114,73,119,78,83,95,49,49,99,104,97,114,95,116,114,97,105,116,115,73,119,69,69,69,69,69,69,0,0,0,0,78,83,116,51,95,95,49,56,116,105,109,101,95,103,101,116,73,99,78,83,95,49,57,105,115,116,114,101,97,109,98,117,102,95,105,116,101,114,97,116,111,114,73,99,78,83,95,49,49,99,104,97,114,95,116,114,97,105,116,115,73,99,69,69,69,69,69,69,0,0,0,0,78,83,116,51,95,95,49,56,110,117,109,112,117,110,99,116,73,119,69,69,0,0,0,0,78,83,116,51,95,95,49,56,110,117,109,112,117,110,99,116,73,99,69,69,0,0,0,0,78,83,116,51,95,95,49,56,109,101,115,115,97,103,101,115,73,119,69,69,0,0,0,0,78,83,116,51,95,95,49,56,109,101,115,115,97,103,101,115,73,99,69,69,0,0,0,0,78,83,116,51,95,95,49,56,105,111,115,95,98,97,115,101,69,0,0,0,0,0,0,0,78,83,116,51,95,95,49,56,105,111,115,95,98,97,115,101,55,102,97,105,108,117,114,101,69,0,0,0,0,0,0,0,78,83,116,51,95,95,49,55,110,117,109,95,112,117,116,73,119,78,83,95,49,57,111,115,116,114,101,97,109,98,117,102,95,105,116,101,114,97,116,111,114,73,119,78,83,95,49,49,99,104,97,114,95,116,114,97,105,116,115,73,119,69,69,69,69,69,69,0,0,0,0,0,78,83,116,51,95,95,49,55,110,117,109,95,112,117,116,73,99,78,83,95,49,57,111,115,116,114,101,97,109,98,117,102,95,105,116,101,114,97,116,111,114,73,99,78,83,95,49,49,99,104,97,114,95,116,114,97,105,116,115,73,99,69,69,69,69,69,69,0,0,0,0,0,78,83,116,51,95,95,49,55,110,117,109,95,103,101,116,73,119,78,83,95,49,57,105,115,116,114,101,97,109,98,117,102,95,105,116,101,114,97,116,111,114,73,119,78,83,95,49,49,99,104,97,114,95,116,114,97,105,116,115,73,119,69,69,69,69,69,69,0,0,0,0,0,78,83,116,51,95,95,49,55,110,117,109,95,103,101,116,73,99,78,83,95,49,57,105,115,116,114,101,97,109,98,117,102,95,105,116,101,114,97,116,111,114,73,99,78,83,95,49,49,99,104,97,114,95,116,114,97,105,116,115,73,99,69,69,69,69,69,69,0,0,0,0,0,78,83,116,51,95,95,49,55,99,111,108,108,97,116,101,73,119,69,69,0,0,0,0,0,78,83,116,51,95,95,49,55,99,111,108,108,97,116,101,73,99,69,69,0,0,0,0,0,78,83,116,51,95,95,49,55,99,111,100,101,99,118,116,73,119,99,49,48,95,109,98,115,116,97,116,101,95,116,69,69,0,0,0,0,0,0,0,0,78,83,116,51,95,95,49,55,99,111,100,101,99,118,116,73,99,99,49,48,95,109,98,115,116,97,116,101,95,116,69,69,0,0,0,0,0,0,0,0,78,83,116,51,95,95,49,55,99,111,100,101,99,118,116,73,68,115,99,49,48,95,109,98,115,116,97,116,101,95,116,69,69,0,0,0,0,0,0,0,78,83,116,51,95,95,49,55,99,111,100,101,99,118,116,73,68,105,99,49,48,95,109,98,115,116,97,116,101,95,116,69,69,0,0,0,0,0,0,0,78,83,116,51,95,95,49,54,108,111,99,97,108,101,53,102,97,99,101,116,69,0,0,0,78,83,116,51,95,95,49,54,108,111,99,97,108,101,53,95,95,105,109,112,69,0,0,0,78,83,116,51,95,95,49,53,99,116,121,112,101,73,119,69,69,0,0,0,0,0,0,0,78,83,116,51,95,95,49,53,99,116,121,112,101,73,99,69,69,0,0,0,0,0,0,0,78,83,116,51,95,95,49,50,49,95,95,98,97,115,105,99,95,115,116,114,105,110,103,95,99,111,109,109,111,110,73,76,98,49,69,69,69,0,0,0,78,83,116,51,95,95,49,50,48,95,95,116,105,109,101,95,103,101,116,95,99,95,115,116,111,114,97,103,101,73,119,69,69,0,0,0,0,0,0,0,78,83,116,51,95,95,49,50,48,95,95,116,105,109,101,95,103,101,116,95,99,95,115,116,111,114,97,103,101,73,99,69,69,0,0,0,0,0,0,0,78,83,116,51,95,95,49,49,57,95,95,105,111,115,116,114,101,97,109,95,99,97,116,101,103,111,114,121,69,0,0,0,78,83,116,51,95,95,49,49,55,95,95,119,105,100,101,110,95,102,114,111,109,95,117,116,102,56,73,76,106,51,50,69,69,69,0,0,0,0,0,0,78,83,116,51,95,95,49,49,54,95,95,110,97,114,114,111,119,95,116,111,95,117,116,102,56,73,76,106,51,50,69,69,69,0,0,0,0,0,0,0,78,83,116,51,95,95,49,49,53,98,97,115,105,99,95,115,116,114,101,97,109,98,117,102,73,119,78,83,95,49,49,99,104,97,114,95,116,114,97,105,116,115,73,119,69,69,69,69,0,0,0,0,0,0,0,0,78,83,116,51,95,95,49,49,53,98,97,115,105,99,95,115,116,114,101,97,109,98,117,102,73,99,78,83,95,49,49,99,104,97,114,95,116,114,97,105,116,115,73,99,69,69,69,69,0,0,0,0,0,0,0,0,78,83,116,51,95,95,49,49,52,101,114,114,111,114,95,99,97,116,101,103,111,114,121,69,0,0,0,0,0,0,0,0,78,83,116,51,95,95,49,49,52,95,95,115,104,97,114,101,100,95,99,111,117,110,116,69,0,0,0,0,0,0,0,0,78,83,116,51,95,95,49,49,52,95,95,110,117,109,95,112,117,116,95,98,97,115,101,69,0,0,0,0,0,0,0,0,78,83,116,51,95,95,49,49,52,95,95,110,117,109,95,103,101,116,95,98,97,115,101,69,0,0,0,0,0,0,0,0,78,83,116,51,95,95,49,49,51,109,101,115,115,97,103,101,115,95,98,97,115,101,69,0,78,83,116,51,95,95,49,49,51,98,97,115,105,99,95,111,115,116,114,101,97,109,73,119,78,83,95,49,49,99,104,97,114,95,116,114,97,105,116,115,73,119,69,69,69,69,0,0,78,83,116,51,95,95,49,49,51,98,97,115,105,99,95,111,115,116,114,101,97,109,73,99,78,83,95,49,49,99,104,97,114,95,116,114,97,105,116,115,73,99,69,69,69,69,0,0,78,83,116,51,95,95,49,49,51,98,97,115,105,99,95,105,115,116,114,101,97,109,73,119,78,83,95,49,49,99,104,97,114,95,116,114,97,105,116,115,73,119,69,69,69,69,0,0,78,83,116,51,95,95,49,49,51,98,97,115,105,99,95,105,115,116,114,101,97,109,73,99,78,83,95,49,49,99,104,97,114,95,116,114,97,105,116,115,73,99,69,69,69,69,0,0,78,83,116,51,95,95,49,49,50,115,121,115,116,101,109,95,101,114,114,111,114,69,0,0,78,83,116,51,95,95,49,49,50,99,111,100,101,99,118,116,95,98,97,115,101,69,0,0,78,83,116,51,95,95,49,49,50,98,97,115,105,99,95,115,116,114,105,110,103,73,119,78,83,95,49,49,99,104,97,114,95,116,114,97,105,116,115,73,119,69,69,78,83,95,57,97,108,108,111,99,97,116,111,114,73,119,69,69,69,69,0,0,78,83,116,51,95,95,49,49,50,98,97,115,105,99,95,115,116,114,105,110,103,73,99,78,83,95,49,49,99,104,97,114,95,116,114,97,105,116,115,73,99,69,69,78,83,95,57,97,108,108,111,99,97,116,111,114,73,99,69,69,69,69,0,0,78,83,116,51,95,95,49,49,50,95,95,100,111,95,109,101,115,115,97,103,101,69,0,0,78,83,116,51,95,95,49,49,49,95,95,115,116,100,111,117,116,98,117,102,73,119,69,69,0,0,0,0,0,0,0,0,78,83,116,51,95,95,49,49,49,95,95,115,116,100,111,117,116,98,117,102,73,99,69,69,0,0,0,0,0,0,0,0,78,83,116,51,95,95,49,49,49,95,95,109,111,110,101,121,95,112,117,116,73,119,69,69,0,0,0,0,0,0,0,0,78,83,116,51,95,95,49,49,49,95,95,109,111,110,101,121,95,112,117,116,73,99,69,69,0,0,0,0,0,0,0,0,78,83,116,51,95,95,49,49,49,95,95,109,111,110,101,121,95,103,101,116,73,119,69,69,0,0,0,0,0,0,0,0,78,83,116,51,95,95,49,49,49,95,95,109,111,110,101,121,95,103,101,116,73,99,69,69,0,0,0,0,0,0,0,0,78,83,116,51,95,95,49,49,48,109,111,110,101,121,112,117,110,99,116,73,119,76,98,49,69,69,69,0,0,0,0,0,78,83,116,51,95,95,49,49,48,109,111,110,101,121,112,117,110,99,116,73,119,76,98,48,69,69,69,0,0,0,0,0,78,83,116,51,95,95,49,49,48,109,111,110,101,121,112,117,110,99,116,73,99,76,98,49,69,69,69,0,0,0,0,0,78,83,116,51,95,95,49,49,48,109,111,110,101,121,112,117,110,99,116,73,99,76,98,48,69,69,69,0,0,0,0,0,78,83,116,51,95,95,49,49,48,109,111,110,101,121,95,98,97,115,101,69,0,0,0,0,78,83,116,51,95,95,49,49,48,99,116,121,112,101,95,98,97,115,101,69,0,0,0,0,78,83,116,51,95,95,49,49,48,95,95,116,105,109,101,95,112,117,116,69,0,0,0,0,78,83,116,51,95,95,49,49,48,95,95,115,116,100,105,110,98,117,102,73,119,69,69,0,78,83,116,51,95,95,49,49,48,95,95,115,116,100,105,110,98,117,102,73,99,69,69,0,78,49,48,101,109,115,99,114,105,112,116,101,110,51,118,97,108,69,0,0,0,0,0,0,78,49,48,101,109,115,99,114,105,112,116,101,110,49,49,109,101,109,111,114,121,95,118,105,101,119,69,0,0,0,0,0,78,49,48,95,95,99,120,120,97,98,105,118,49,50,51,95,95,102,117,110,100,97,109,101,110,116,97,108,95,116,121,112,101,95,105,110,102,111,69,0,78,49,48,95,95,99,120,120,97,98,105,118,49,50,49,95,95,118,109,105,95,99,108,97,115,115,95,116,121,112,101,95,105,110,102,111,69,0,0,0,78,49,48,95,95,99,120,120,97,98,105,118,49,50,48,95,95,115,105,95,99,108,97,115,115,95,116,121,112,101,95,105,110,102,111,69,0,0,0,0,78,49,48,95,95,99,120,120,97,98,105,118,49,49,57,95,95,112,111,105,110,116,101,114,95,116,121,112,101,95,105,110,102,111,69,0,0,0,0,0,78,49,48,95,95,99,120,120,97,98,105,118,49,49,55,95,95,112,98,97,115,101,95,116,121,112,101,95,105,110,102,111,69,0,0,0,0,0,0,0,78,49,48,95,95,99,120,120,97,98,105,118,49,49,55,95,95,99,108,97,115,115,95,116,121,112,101,95,105,110,102,111,69,0,0,0,0,0,0,0,78,49,48,95,95,99,120,120,97,98,105,118,49,49,54,95,95,115,104,105,109,95,116,121,112,101,95,105,110,102,111,69,0,0,0,0,0,0,0,0,68,110,0,0,0,0,0,0,56,65,103,114,97,112,104,95,115,0,0,0,0,0,0,0,56,65,103,110,111,100,101,95,115,0,0,0,0,0,0,0,56,65,103,101,100,103,101,95,115,0,0,0,0,0,0,0,56,65,103,100,101,115,99,95,115,0,0,0,0,0,0,0,53,71,86,67,95,115,0,0,200,203,2,0,40,204,2,0,200,203,2,0,136,204,2,0,0,0,0,0,152,204,2,0,0,0,0,0,168,204,2,0,0,0,0,0,184,204,2,0,160,218,2,0,0,0,0,0,0,0,0,0,200,204,2,0,160,218,2,0,0,0,0,0,0,0,0,0,216,204,2,0,160,218,2,0,0,0,0,0,0,0,0,0,240,204,2,0,248,218,2,0,0,0,0,0,0,0,0,0,8,205,2,0,248,218,2,0,0,0,0,0,0,0,0,0,32,205,2,0,160,218,2,0,0,0,0,0,0,0,0,0,56,205,2,0,0,0,0,0,0,0,0,0,0,0,0,0,64,205,2,0,1,0,0,0,224,225,2,0,0,0,0,0,80,205,2,0,1,0,0,0,232,225,2,0,0,0,0,0,96,205,2,0,1,0,0,0,240,225,2,0,0,0,0,0,112,205,2,0,1,0,0,0,0,226,2,0,0,0,0,0,128,205,2,0,0,0,0,0,224,225,2,0,0,0,0,0,144,205,2,0,0,0,0,0,232,225,2,0,0,0,0,0,160,205,2,0,0,0,0,0,240,225,2,0,0,0,0,0,176,205,2,0,0,0,0,0,0,226,2,0,0,0,0,0,184,205,2,0,240,203,2,0,208,205,2,0,0,0,0,0,2,0,0,0,200,222,2,0,2,0,0,0,128,224,2,0,0,0,0,0,240,203,2,0,24,206,2,0,0,0,0,0,2,0,0,0,200,222,2,0,2,0,0,0,136,224,2,0,0,0,0,0,240,203,2,0,96,206,2,0,0,0,0,0,2,0,0,0,200,222,2,0,2,0,0,0,144,224,2,0,0,0,0,0,240,203,2,0,168,206,2,0,0,0,0,0,2,0,0,0,200,222,2,0,2,0,0,0,152,224,2,0,0,0,0,0,0,0,0,0,240,206,2,0,144,221,2,0,0,0,0,0,0,0,0,0,32,207,2,0,144,221,2,0,0,0,0,0,240,203,2,0,80,207,2,0,0,0,0,0,1,0,0,0,144,223,2,0,0,0,0,0,240,203,2,0,104,207,2,0,0,0,0,0,1,0,0,0,144,223,2,0,0,0,0,0,240,203,2,0,128,207,2,0,0,0,0,0,1,0,0,0,152,223,2,0,0,0,0,0,240,203,2,0,152,207,2,0,0,0,0,0,1,0,0,0,152,223,2,0,0,0,0,0,240,203,2,0,176,207,2,0,0,0,0,0,2,0,0,0,200,222,2,0,2,0,0,0,48,225,2,0,0,8,0,0,240,203,2,0,248,207,2,0,0,0,0,0,2,0,0,0,200,222,2,0,2,0,0,0,48,225,2,0,0,8,0,0,240,203,2,0,64,208,2,0,0,0,0,0,3,0,0,0,200,222,2,0,2,0,0,0,152,219,2,0,2,0,0,0,48,223,2,0,0,8,0,0,240,203,2,0,136,208,2,0,0,0,0,0,3,0,0,0,200,222,2,0,2,0,0,0,152,219,2,0,2,0,0,0,56,223,2,0,0,8,0,0,0,0,0,0,208,208,2,0,200,222,2,0,0,0,0,0,0,0,0,0,232,208,2,0,200,222,2,0,0,0,0,0,240,203,2,0,0,209,2,0,0,0,0,0,2,0,0,0,200,222,2,0,2,0,0,0,160,223,2,0,2,0,0,0,240,203,2,0,24,209,2,0,0,0,0,0,2,0,0,0,200,222,2,0,2,0,0,0,160,223,2,0,2,0,0,0,0,0,0,0,48,209,2,0,0,0,0,0,72,209,2,0,8,224,2,0,0,0,0,0,240,203,2,0,104,209,2,0,0,0,0,0,2,0,0,0,200,222,2,0,2,0,0,0,64,220,2,0,0,0,0,0,240,203,2,0,176,209,2,0,0,0,0,0,2,0,0,0,200,222,2,0,2,0,0,0,88,220,2,0,0,0,0,0,240,203,2,0,248,209,2,0,0,0,0,0,2,0,0,0,200,222,2,0,2,0,0,0,112,220,2,0,0,0,0,0,240,203,2,0,64,210,2,0,0,0,0,0,2,0,0,0,200,222,2,0,2,0,0,0,136,220,2,0,0,0,0,0,0,0,0,0,136,210,2,0,200,222,2,0,0,0,0,0,0,0,0,0,160,210,2,0,200,222,2,0,0,0,0,0,240,203,2,0,184,210,2,0,0,0,0,0,2,0,0,0,200,222,2,0,2,0,0,0,24,224,2,0,2,0,0,0,240,203,2,0,224,210,2,0,0,0,0,0,2,0,0,0,200,222,2,0,2,0,0,0,24,224,2,0,2,0,0,0,240,203,2,0,8,211,2,0,0,0,0,0,2,0,0,0,200,222,2,0,2,0,0,0,24,224,2,0,2,0,0,0,240,203,2,0,48,211,2,0,0,0,0,0,2,0,0,0,200,222,2,0,2,0,0,0,24,224,2,0,2,0,0,0,0,0,0,0,88,211,2,0,136,223,2,0,0,0,0,0,0,0,0,0,112,211,2,0,200,222,2,0,0,0,0,0,240,203,2,0,136,211,2,0,0,0,0,0,2,0,0,0,200,222,2,0,2,0,0,0,40,225,2,0,2,0,0,0,240,203,2,0,160,211,2,0,0,0,0,0,2,0,0,0,200,222,2,0,2,0,0,0,40,225,2,0,2,0,0,0,0,0,0,0,184,211,2,0,0,0,0,0,224,211,2,0,0,0,0,0,8,212,2,0,0,0,0,0,48,212,2,0,80,224,2,0,0,0,0,0,0,0,0,0,80,212,2,0,168,222,2,0,0,0,0,0,0,0,0,0,120,212,2,0,168,222,2,0,0,0,0,0,0,0,0,0,160,212,2,0,0,0,0,0,216,212,2,0,0,0,0,0,16,213,2,0,0,0,0,0,48,213,2,0,0,0,0,0,80,213,2,0,0,0,0,0,112,213,2,0,0,0,0,0,144,213,2,0,240,203,2,0,168,213,2,0,0,0,0,0,1,0,0,0,32,220,2,0,3,244,255,255,240,203,2,0,216,213,2,0,0,0,0,0,1,0,0,0,48,220,2,0,3,244,255,255,240,203,2,0,8,214,2,0,0,0,0,0,1,0,0,0,32,220,2,0,3,244,255,255,240,203,2,0,56,214,2,0,0,0,0,0,1,0,0,0,48,220,2,0,3,244,255,255,0,0,0,0,104,214,2,0,200,218,2,0,0,0,0,0,0,0,0,0,128,214,2,0,240,203,2,0,152,214,2,0,0,0,0,0,1,0,0,0,40,223,2,0,0,0,0,0,240,203,2,0,216,214,2,0,0,0,0,0,1,0,0,0,40,223,2,0,0,0,0,0,0,0,0,0,24,215,2,0,128,223,2,0,0,0,0,0,0,0,0,0,48,215,2,0,112,223,2,0,0,0,0,0,0,0,0,0,80,215,2,0,120,223,2,0,0,0,0,0,0,0,0,0,112,215,2,0,0,0,0,0,144,215,2,0,0,0,0,0,176,215,2,0,0,0,0,0,208,215,2,0,240,203,2,0,240,215,2,0,0,0,0,0,2,0,0,0,200,222,2,0,2,0,0,0,32,225,2,0,2,0,0,0,240,203,2,0,16,216,2,0,0,0,0,0,2,0,0,0,200,222,2,0,2,0,0,0,32,225,2,0,2,0,0,0,240,203,2,0,48,216,2,0,0,0,0,0,2,0,0,0,200,222,2,0,2,0,0,0,32,225,2,0,2,0,0,0,240,203,2,0,80,216,2,0,0,0,0,0,2,0,0,0,200,222,2,0,2,0,0,0,32,225,2,0,2,0,0,0,0,0,0,0,112,216,2,0,0,0,0,0,136,216,2,0,0,0,0,0,160,216,2,0,0,0,0,0,184,216,2,0,112,223,2,0,0,0,0,0,0,0,0,0,208,216,2,0,120,223,2,0,0,0,0,0,0,0,0,0,232,216,2,0,0,0,0,0,0,217,2,0,0,0,0,0,32,217,2,0,200,225,2,0,0,0,0,0,0,0,0,0,72,217,2,0,184,225,2,0,0,0,0,0,0,0,0,0,112,217,2,0,184,225,2,0,0,0,0,0,0,0,0,0,152,217,2,0,168,225,2,0,0,0,0,0,0,0,0,0,192,217,2,0,200,225,2,0,0,0,0,0,0,0,0,0,232,217,2,0,200,225,2,0,0,0,0,0,0,0,0,0,16,218,2,0,152,218,2,0,0,0,0,0,200,203,2,0,56,218,2,0,0,0,0,0,64,218,2,0,0,0,0,0,80,218,2,0,0,0,0,0,96,218,2,0,0,0,0,0,112,218,2,0,0,0,0,0,128,218,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,255,255,255,255,0,0,0,0,255,255,255,255,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,48,49,50,51,52,53,54,55,56,57,97,98,99,100,101,102,65,66,67,68,69,70,120,88,43,45,112,80,105,73,110,78,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,30,4,0,0,4,0,0,0,208,2,0,0,64,0,0,0,30,4,0,0,4,0,0,0,30,4,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,168,5,2,0,16,33,0,0,184,44,0,0,0,0,0,0,200,108,2,0,16,33,0,0,216,39,0,0,0,0,0,0,8,154,2,0,16,33,0,0,8,43,0,0,0,0,0,0,96,110,2,0,16,33,0,0,8,43,0,0,0,0,0,0,184,97,2,0,16,33,0,0,40,44,0,0,0,0,0,0,224,84,2,0,88,33,0,0,40,44,0,0,0,0,0,0,216,73,2,0,16,33,0,0,56,43,0,0,0,0,0,0,48,57,2,0,16,33,0,0,168,36,0,0,0,0,0,0,136,218,1,0,16,33,0,0,8,40,0,0,0,0,0,0,80,41,2,0,16,33,0,0,8,40,0,0,0,0,0,0,0,35,2,0,16,33,0,0,200,43,0,0,0,0,0,0,96,28,2,0,16,33,0,0,216,36,0,0,0,0,0,0,248,19,2,0,16,33,0,0,104,40,0,0,0,0,0,0,96,10,2,0,16,33,0,0,72,42,0,0,0,0,0,0,176,3,2,0,16,33,0,0,56,40,0,0,0,0,0,0,56,1,2,0,16,33,0,0,120,42,0,0,0,0,0,0,200,254,1,0,16,33,0,0,248,37,0,0,0,0,0,0,232,252,1,0,16,33,0,0,152,40,0,0,0,0,0,0,0,251,1,0,16,33,0,0,248,40,0,0,0,0,0,0,96,159,2,0,16,33,0,0,104,37,0,0,0,0,0,0,16,246,1,0,16,33,0,0,168,42,0,0,0,0,0,0,240,243,1,0,16,33,0,0,136,44,0,0,0,0,0,0,208,240,1,0,16,33,0,0,248,43,0,0,0,0,0,0,152,237,1,0,16,33,0,0,184,44,0,0,0,0,0,0,168,234,1,0,16,33,0,0,184,44,0,0,0,0,0,0,80,232,1,0,16,33,0,0,152,37,0,0,0,0,0,0,48,230,1,0,16,33,0,0,152,43,0,0,0,0,0,0,168,227,1,0,16,33,0,0,104,43,0,0,0,0,0,0,104,225,1,0,16,33,0,0,120,36,0,0,0,0,0,0,208,222,1,0,16,33,0,0,136,41,0,0,0,0,0,0,48,221,1,0,16,33,0,0,184,41,0,0,0,0,0,0,72,219,1,0,16,33,0,0,232,41,0,0,0,0,0,0,56,217,1,0,16,33,0,0,72,45,0,0,0,0,0,0,80,213,1,0,16,33,0,0,24,45,0,0,0,0,0,0,200,210,1,0,16,33,0,0,120,45,0,0,0,0,0,0,128,208,1,0,16,33,0,0,120,39,0,0,0,0,0,0,120,206,1,0,16,33,0,0,88,44,0,0,0,0,0,0,176,204,1,0,16,33,0,0,56,37,0,0,0,0,0,0,240,202,1,0,16,33,0,0,72,36,0,0,0,0,0,0,32,201,1,0,16,33,0,0,24,42,0,0,0,0,0,0,184,199,1,0,16,33,0,0,136,38,0,0,0,0,0,0,0,198,1,0,16,33,0,0,88,38,0,0,0,0,0,0,8,195,1,0,16,33,0,0,72,39,0,0,0,0,0,0,32,192,1,0,16,33,0,0,24,39,0,0,0,0,0,0,64,190,1,0,16,33,0,0,168,39,0,0,0,0,0,0,0,188,1,0,16,33,0,0,184,38,0,0,0,0,0,0,152,186,1,0,16,33,0,0,216,42,0,0,0,0,0,0,240,184,1,0,16,33,0,0,8,37,0,0,0,0,0,0,56,183,1,0,16,33,0,0,200,40,0,0,0,0,0,0,136,181,1,0,16,33,0,0,232,44,0,0,0,0,0,0,192,179,1,0,16,33,0,0,200,37,0,0,0,0,0,0,208,177,1,0,16,33,0,0,40,38,0,0,0,0,0,0,160,175,1,0,16,33,0,0,88,41,0,0,0,0,0,0,224,172,1,0,16,33,0,0,232,38,0,0,0,0,0,0,232,169,1,0,16,33,0,0,40,41,0,0,0,0,0,0,56,168,1,0,224,16,0,0,0,0,0,0,0,0,0,0,32,166,1,0,224,16,0,0,0,0,0,0,0,0,0,0,160,217,1,0,136,74,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,12,0,0,0,255,255,255,255,0,0,0,0,0,0,0,0,24,0,0,0,0,0,0,0,0,0,0,0,126,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,53,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,232,4,0,0,0,0,0,0,126,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,250,3,0,0,0,0,0,0,126,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,232,230,1,0,88,119,2,0,224,15,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,107,101,121,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,0,0,0,255,255,255,255,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,0,0,0,255,255,255,255,0,0,0,0,0,0,0,0,68,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,0,0,0,4,0,0,0,0,0,0,0,102,3,0,0,236,0,0,0,248,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,204,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,95,65,71,95,100,97,116,97,100,105,99,116,0,0,0,0,0,0,0,0,0,0,0,0,95,65,71,95,112,101,110,100,105,110,103,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,240,191,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,240,63,250,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,240,63,92,4,0,0,0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,224,63,214,4,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,240,63,8,0,0,0,0,0,0,0,5,0,0,0,0,0,0,0,51,51,51,51,51,51,243,63,148,0,0,0,0,0,0,0,6,0,0,0,0,0,0,0,154,153,153,153,153,153,233,63,46,0,0,0,0,0,0,0,7,0,0,0,0,0,0,0,0,0,0,0,0,0,240,63,182,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,80,65,2,0,25,0,0,0,0,0,0,0,0,0,0,0,152,122,2,0,1,0,0,0,224,87,2,0,2,0,0,0,152,1,2,0,3,0,0,0,168,5,2,0,4,0,0,0,0,35,2,0,5,0,0,0,80,235,1,0,6,0,0,0,136,218,1,0,0,0,0,0,88,148,1,0,17,0,0,0,120,131,1,0,18,0,0,0,224,114,1,0,18,0,0,0,216,163,2,0,1,0,0,0,40,149,2,0,7,0,0,0,0,0,0,0,0,0,0,0,40,128,1,0,8,0,0,0,104,113,1,0,64,0,0,0,72,204,2,0,32,0,0,0,112,204,2,0,8,0,0,0,136,80,2,0,32,0,0,0,0,0,0,0,0,0,0,0,192,140,1,0,0,0,0,0,1,0,0,0,88,139,1,0,1,0,0,0,0,0,0,0,232,136,1,0,1,0,0,0,1,0,0,0,136,218,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,0,0,0,0,0,0,0,10,0,0,0,0,0,0,0,11,0,0,0,0,0,0,0,9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,118,4,0,0,178,1,0,0,0,0,0,0,126,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,58,1,0,0,0,0,0,0,126,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,16,0,0,0,0,0,0,0,0,0,0,0,96,3,0,0,0,0,0,0,126,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,255,255,255,255,0,0,0,0,0,0,0,0,208,0,0,0,0,0,0,0,126,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,255,255,255,255,0,0,0,0,0,0,0,0,106,2,0,0,0,0,0,0,126,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,20,0,0,0,0,0,0,0,0,0,0,0,208,0,0,0,0,0,0,0,126,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,12,0,0,0,0,0,0,0,0,0,0,0,106,2,0,0,0,0,0,0,126,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,184,3,0,0,30,3,0,0,180,0,0,0,20,1,0,0,0,0,0,0,0,0,0,0,138,1,0,0,224,2,0,0,54,5,0,0,0,0,0,0,0,4,0,0,6,5,0,0,78,3,0,0,174,5,0,0,126,4,0,0,26,3,0,0,216,3,0,0,0,0,0,0,112,245,2,0,152,245,2,0,136,245,2,0,0,0,0,0,8,0,0,0,255,255,255,255,0,0,0,0,0,0,0,0,0,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
, "i8", ALLOC_NONE, Runtime.GLOBAL_BASE)
function runPostSets() {
HEAP32[(((__ZTVN10__cxxabiv120__si_class_type_infoE)+(8))>>2)]=(1422);
HEAP32[(((__ZTVN10__cxxabiv120__si_class_type_infoE)+(12))>>2)]=(700);
HEAP32[(((__ZTVN10__cxxabiv120__si_class_type_infoE)+(16))>>2)]=(528);
HEAP32[(((__ZTVN10__cxxabiv120__si_class_type_infoE)+(20))>>2)]=(520);
HEAP32[(((__ZTVN10__cxxabiv120__si_class_type_infoE)+(24))>>2)]=(430);
HEAP32[(((__ZTVN10__cxxabiv120__si_class_type_infoE)+(28))>>2)]=(210);
HEAP32[(((__ZTVN10__cxxabiv120__si_class_type_infoE)+(32))>>2)]=(498);
HEAP32[(((__ZTVN10__cxxabiv120__si_class_type_infoE)+(36))>>2)]=(548);
HEAP32[(((__ZTVN10__cxxabiv119__pointer_type_infoE)+(8))>>2)]=(1422);
HEAP32[(((__ZTVN10__cxxabiv119__pointer_type_infoE)+(12))>>2)]=(354);
HEAP32[(((__ZTVN10__cxxabiv119__pointer_type_infoE)+(16))>>2)]=(528);
HEAP32[(((__ZTVN10__cxxabiv119__pointer_type_infoE)+(20))>>2)]=(520);
HEAP32[(((__ZTVN10__cxxabiv119__pointer_type_infoE)+(24))>>2)]=(990);
HEAP32[(((__ZTVN10__cxxabiv117__class_type_infoE)+(8))>>2)]=(1422);
HEAP32[(((__ZTVN10__cxxabiv117__class_type_infoE)+(12))>>2)]=(1408);
HEAP32[(((__ZTVN10__cxxabiv117__class_type_infoE)+(16))>>2)]=(528);
HEAP32[(((__ZTVN10__cxxabiv117__class_type_infoE)+(20))>>2)]=(520);
HEAP32[(((__ZTVN10__cxxabiv117__class_type_infoE)+(24))>>2)]=(430);
HEAP32[(((__ZTVN10__cxxabiv117__class_type_infoE)+(28))>>2)]=(1100);
HEAP32[(((__ZTVN10__cxxabiv117__class_type_infoE)+(32))>>2)]=(542);
HEAP32[(((__ZTVN10__cxxabiv117__class_type_infoE)+(36))>>2)]=(830);
HEAP32[((__ZTIt)>>2)]=(((183240)|0));
HEAP32[(((__ZTIt)+(4))>>2)]=((183344)|0);
HEAP32[((__ZTIs)>>2)]=(((183240)|0));
HEAP32[(((__ZTIs)+(4))>>2)]=((183352)|0);
HEAP32[((__ZTIm)>>2)]=(((183240)|0));
HEAP32[(((__ZTIm)+(4))>>2)]=((183360)|0);
HEAP32[((__ZTIl)>>2)]=(((183240)|0));
HEAP32[(((__ZTIl)+(4))>>2)]=((183368)|0);
HEAP32[((__ZTIj)>>2)]=(((183240)|0));
HEAP32[(((__ZTIj)+(4))>>2)]=((183376)|0);
HEAP32[((__ZTIi)>>2)]=(((183240)|0));
HEAP32[(((__ZTIi)+(4))>>2)]=((183384)|0);
HEAP32[((__ZTIh)>>2)]=(((183240)|0));
HEAP32[(((__ZTIh)+(4))>>2)]=((183392)|0);
HEAP32[((__ZTIf)>>2)]=(((183240)|0));
HEAP32[(((__ZTIf)+(4))>>2)]=((183400)|0);
HEAP32[((__ZTId)>>2)]=(((183240)|0));
HEAP32[(((__ZTId)+(4))>>2)]=((183416)|0);
HEAP32[((__ZTIc)>>2)]=(((183240)|0));
HEAP32[(((__ZTIc)+(4))>>2)]=((183424)|0);
HEAP32[((__ZTIa)>>2)]=(((183240)|0));
HEAP32[(((__ZTIa)+(4))>>2)]=((183440)|0);
HEAP32[((187032)>>2)]=(((__ZTVN10__cxxabiv117__class_type_infoE+8)|0));
HEAP32[((187040)>>2)]=(((__ZTVN10__cxxabiv117__class_type_infoE+8)|0));
HEAP32[((187048)>>2)]=(((__ZTVN10__cxxabiv120__si_class_type_infoE+8)|0));
HEAP32[((187064)>>2)]=(((__ZTVN10__cxxabiv120__si_class_type_infoE+8)|0));
HEAP32[((187080)>>2)]=(((__ZTVN10__cxxabiv120__si_class_type_infoE+8)|0));
HEAP32[((187096)>>2)]=(((__ZTVN10__cxxabiv120__si_class_type_infoE+8)|0));
HEAP32[((187112)>>2)]=(((__ZTVN10__cxxabiv120__si_class_type_infoE+8)|0));
HEAP32[((187128)>>2)]=(((__ZTVN10__cxxabiv120__si_class_type_infoE+8)|0));
HEAP32[((187144)>>2)]=(((__ZTVN10__cxxabiv119__pointer_type_infoE+8)|0));
HEAP32[((187156)>>2)]=__ZTIc;
HEAP32[((187160)>>2)]=(((__ZTVN10__cxxabiv119__pointer_type_infoE+8)|0));
HEAP32[((187176)>>2)]=(((__ZTVN10__cxxabiv119__pointer_type_infoE+8)|0));
HEAP32[((187192)>>2)]=(((__ZTVN10__cxxabiv119__pointer_type_infoE+8)|0));
HEAP32[((187208)>>2)]=(((__ZTVN10__cxxabiv119__pointer_type_infoE+8)|0));
HEAP32[((187224)>>2)]=(((__ZTVN10__cxxabiv119__pointer_type_infoE+8)|0));
HEAP32[((187240)>>2)]=(((__ZTVN10__cxxabiv119__pointer_type_infoE+8)|0));
HEAP32[((187256)>>2)]=(((__ZTVN10__cxxabiv119__pointer_type_infoE+8)|0));
HEAP32[((187272)>>2)]=(((__ZTVN10__cxxabiv119__pointer_type_infoE+8)|0));
HEAP32[((187288)>>2)]=(((__ZTVN10__cxxabiv117__class_type_infoE+8)|0));
HEAP32[((187424)>>2)]=(((__ZTVN10__cxxabiv120__si_class_type_infoE+8)|0));
HEAP32[((187440)>>2)]=(((__ZTVN10__cxxabiv120__si_class_type_infoE+8)|0));
HEAP32[((187696)>>2)]=(((__ZTVN10__cxxabiv120__si_class_type_infoE+8)|0));
HEAP32[((187712)>>2)]=(((__ZTVN10__cxxabiv120__si_class_type_infoE+8)|0));
HEAP32[((187792)>>2)]=(((__ZTVN10__cxxabiv117__class_type_infoE+8)|0));
HEAP32[((187800)>>2)]=(((__ZTVN10__cxxabiv120__si_class_type_infoE+8)|0));
HEAP32[((187944)>>2)]=(((__ZTVN10__cxxabiv120__si_class_type_infoE+8)|0));
HEAP32[((187960)>>2)]=(((__ZTVN10__cxxabiv120__si_class_type_infoE+8)|0));
HEAP32[((188104)>>2)]=(((__ZTVN10__cxxabiv120__si_class_type_infoE+8)|0));
HEAP32[((188120)>>2)]=(((__ZTVN10__cxxabiv120__si_class_type_infoE+8)|0));
HEAP32[((188200)>>2)]=(((__ZTVN10__cxxabiv117__class_type_infoE+8)|0));
HEAP32[((188208)>>2)]=(((__ZTVN10__cxxabiv117__class_type_infoE+8)|0));
HEAP32[((188216)>>2)]=(((__ZTVN10__cxxabiv117__class_type_infoE+8)|0));
HEAP32[((188224)>>2)]=(((__ZTVN10__cxxabiv120__si_class_type_infoE+8)|0));
HEAP32[((188240)>>2)]=(((__ZTVN10__cxxabiv120__si_class_type_infoE+8)|0));
HEAP32[((188256)>>2)]=(((__ZTVN10__cxxabiv120__si_class_type_infoE+8)|0));
HEAP32[((188272)>>2)]=(((__ZTVN10__cxxabiv117__class_type_infoE+8)|0));
HEAP32[((188280)>>2)]=(((__ZTVN10__cxxabiv117__class_type_infoE+8)|0));
HEAP32[((188288)>>2)]=(((__ZTVN10__cxxabiv117__class_type_infoE+8)|0));
HEAP32[((188296)>>2)]=(((__ZTVN10__cxxabiv117__class_type_infoE+8)|0));
HEAP32[((188304)>>2)]=(((__ZTVN10__cxxabiv117__class_type_infoE+8)|0));
HEAP32[((188312)>>2)]=(((__ZTVN10__cxxabiv117__class_type_infoE+8)|0));
HEAP32[((188320)>>2)]=(((__ZTVN10__cxxabiv117__class_type_infoE+8)|0));
HEAP32[((188424)>>2)]=(((__ZTVN10__cxxabiv120__si_class_type_infoE+8)|0));
HEAP32[((188440)>>2)]=(((__ZTVN10__cxxabiv117__class_type_infoE+8)|0));
HEAP32[((188496)>>2)]=(((__ZTVN10__cxxabiv120__si_class_type_infoE+8)|0));
HEAP32[((188512)>>2)]=(((__ZTVN10__cxxabiv120__si_class_type_infoE+8)|0));
HEAP32[((188528)>>2)]=(((__ZTVN10__cxxabiv120__si_class_type_infoE+8)|0));
HEAP32[((188544)>>2)]=(((__ZTVN10__cxxabiv117__class_type_infoE+8)|0));
HEAP32[((188552)>>2)]=(((__ZTVN10__cxxabiv117__class_type_infoE+8)|0));
HEAP32[((188560)>>2)]=(((__ZTVN10__cxxabiv117__class_type_infoE+8)|0));
HEAP32[((188568)>>2)]=(((__ZTVN10__cxxabiv117__class_type_infoE+8)|0));
HEAP32[((188704)>>2)]=(((__ZTVN10__cxxabiv117__class_type_infoE+8)|0));
HEAP32[((188712)>>2)]=(((__ZTVN10__cxxabiv117__class_type_infoE+8)|0));
HEAP32[((188720)>>2)]=(((__ZTVN10__cxxabiv117__class_type_infoE+8)|0));
HEAP32[((188728)>>2)]=(((__ZTVN10__cxxabiv120__si_class_type_infoE+8)|0));
HEAP32[((188744)>>2)]=(((__ZTVN10__cxxabiv120__si_class_type_infoE+8)|0));
HEAP32[((188760)>>2)]=(((__ZTVN10__cxxabiv117__class_type_infoE+8)|0));
HEAP32[((188768)>>2)]=(((__ZTVN10__cxxabiv117__class_type_infoE+8)|0));
HEAP32[((188776)>>2)]=(((__ZTVN10__cxxabiv120__si_class_type_infoE+8)|0));
HEAP32[((188792)>>2)]=(((__ZTVN10__cxxabiv120__si_class_type_infoE+8)|0));
HEAP32[((188808)>>2)]=(((__ZTVN10__cxxabiv120__si_class_type_infoE+8)|0));
HEAP32[((188824)>>2)]=(((__ZTVN10__cxxabiv120__si_class_type_infoE+8)|0));
HEAP32[((188840)>>2)]=(((__ZTVN10__cxxabiv120__si_class_type_infoE+8)|0));
HEAP32[((188856)>>2)]=(((__ZTVN10__cxxabiv120__si_class_type_infoE+8)|0));
HEAP32[((188872)>>2)]=(((__ZTVN10__cxxabiv120__si_class_type_infoE+8)|0));
HEAP32[((188896)>>2)]=(((__ZTVN10__cxxabiv117__class_type_infoE+8)|0));
HEAP32[((188904)>>2)]=(((__ZTVN10__cxxabiv117__class_type_infoE+8)|0));
HEAP32[((188912)>>2)]=(((__ZTVN10__cxxabiv117__class_type_infoE+8)|0));
HEAP32[((188920)>>2)]=(((__ZTVN10__cxxabiv117__class_type_infoE+8)|0));
HEAP32[((188928)>>2)]=(((__ZTVN10__cxxabiv117__class_type_infoE+8)|0));
}
if (!awaitingMemoryInitializer) runPostSets();
var tempDoublePtr = Runtime.alignMemory(allocate(12, "i8", ALLOC_STATIC), 8);
assert(tempDoublePtr % 8 == 0);
function copyTempFloat(ptr) { // functions, because inlining this code increases code size too much
  HEAP8[tempDoublePtr] = HEAP8[ptr];
  HEAP8[tempDoublePtr+1] = HEAP8[ptr+1];
  HEAP8[tempDoublePtr+2] = HEAP8[ptr+2];
  HEAP8[tempDoublePtr+3] = HEAP8[ptr+3];
}
function copyTempDouble(ptr) {
  HEAP8[tempDoublePtr] = HEAP8[ptr];
  HEAP8[tempDoublePtr+1] = HEAP8[ptr+1];
  HEAP8[tempDoublePtr+2] = HEAP8[ptr+2];
  HEAP8[tempDoublePtr+3] = HEAP8[ptr+3];
  HEAP8[tempDoublePtr+4] = HEAP8[ptr+4];
  HEAP8[tempDoublePtr+5] = HEAP8[ptr+5];
  HEAP8[tempDoublePtr+6] = HEAP8[ptr+6];
  HEAP8[tempDoublePtr+7] = HEAP8[ptr+7];
}
  function ___gxx_personality_v0() {
    }
  function _memcpy(dest, src, num) {
      dest = dest|0; src = src|0; num = num|0;
      var ret = 0;
      ret = dest|0;
      if ((dest&3) == (src&3)) {
        while (dest & 3) {
          if ((num|0) == 0) return ret|0;
          HEAP8[(dest)]=HEAP8[(src)];
          dest = (dest+1)|0;
          src = (src+1)|0;
          num = (num-1)|0;
        }
        while ((num|0) >= 4) {
          HEAP32[((dest)>>2)]=HEAP32[((src)>>2)];
          dest = (dest+4)|0;
          src = (src+4)|0;
          num = (num-4)|0;
        }
      }
      while ((num|0) > 0) {
        HEAP8[(dest)]=HEAP8[(src)];
        dest = (dest+1)|0;
        src = (src+1)|0;
        num = (num-1)|0;
      }
      return ret|0;
    }var _llvm_memcpy_p0i8_p0i8_i32=_memcpy;
  var ERRNO_CODES={EPERM:1,ENOENT:2,ESRCH:3,EINTR:4,EIO:5,ENXIO:6,E2BIG:7,ENOEXEC:8,EBADF:9,ECHILD:10,EAGAIN:11,EWOULDBLOCK:11,ENOMEM:12,EACCES:13,EFAULT:14,ENOTBLK:15,EBUSY:16,EEXIST:17,EXDEV:18,ENODEV:19,ENOTDIR:20,EISDIR:21,EINVAL:22,ENFILE:23,EMFILE:24,ENOTTY:25,ETXTBSY:26,EFBIG:27,ENOSPC:28,ESPIPE:29,EROFS:30,EMLINK:31,EPIPE:32,EDOM:33,ERANGE:34,ENOMSG:35,EIDRM:36,ECHRNG:37,EL2NSYNC:38,EL3HLT:39,EL3RST:40,ELNRNG:41,EUNATCH:42,ENOCSI:43,EL2HLT:44,EDEADLK:45,ENOLCK:46,EBADE:50,EBADR:51,EXFULL:52,ENOANO:53,EBADRQC:54,EBADSLT:55,EDEADLOCK:56,EBFONT:57,ENOSTR:60,ENODATA:61,ETIME:62,ENOSR:63,ENONET:64,ENOPKG:65,EREMOTE:66,ENOLINK:67,EADV:68,ESRMNT:69,ECOMM:70,EPROTO:71,EMULTIHOP:74,ELBIN:75,EDOTDOT:76,EBADMSG:77,EFTYPE:79,ENOTUNIQ:80,EBADFD:81,EREMCHG:82,ELIBACC:83,ELIBBAD:84,ELIBSCN:85,ELIBMAX:86,ELIBEXEC:87,ENOSYS:88,ENMFILE:89,ENOTEMPTY:90,ENAMETOOLONG:91,ELOOP:92,EOPNOTSUPP:95,EPFNOSUPPORT:96,ECONNRESET:104,ENOBUFS:105,EAFNOSUPPORT:106,EPROTOTYPE:107,ENOTSOCK:108,ENOPROTOOPT:109,ESHUTDOWN:110,ECONNREFUSED:111,EADDRINUSE:112,ECONNABORTED:113,ENETUNREACH:114,ENETDOWN:115,ETIMEDOUT:116,EHOSTDOWN:117,EHOSTUNREACH:118,EINPROGRESS:119,EALREADY:120,EDESTADDRREQ:121,EMSGSIZE:122,EPROTONOSUPPORT:123,ESOCKTNOSUPPORT:124,EADDRNOTAVAIL:125,ENETRESET:126,EISCONN:127,ENOTCONN:128,ETOOMANYREFS:129,EPROCLIM:130,EUSERS:131,EDQUOT:132,ESTALE:133,ENOTSUP:134,ENOMEDIUM:135,ENOSHARE:136,ECASECLASH:137,EILSEQ:138,EOVERFLOW:139,ECANCELED:140,ENOTRECOVERABLE:141,EOWNERDEAD:142,ESTRPIPE:143};
  var ___errno_state=0;function ___setErrNo(value) {
      // For convenient setting and returning of errno.
      HEAP32[((___errno_state)>>2)]=value
      return value;
    }
  var _stdin=allocate(1, "i32*", ALLOC_STATIC);
  var _stdout=allocate(1, "i32*", ALLOC_STATIC);
  var _stderr=allocate(1, "i32*", ALLOC_STATIC);
  var __impure_ptr=allocate(1, "i32*", ALLOC_STATIC);var FS={currentPath:"/",nextInode:2,streams:[null],ignorePermissions:true,createFileHandle:function (stream, fd) {
        if (typeof stream === 'undefined') {
          stream = null;
        }
        if (!fd) {
          if (stream && stream.socket) {
            for (var i = 1; i < 64; i++) {
              if (!FS.streams[i]) {
                fd = i;
                break;
              }
            }
            assert(fd, 'ran out of low fds for sockets');
          } else {
            fd = Math.max(FS.streams.length, 64);
            for (var i = FS.streams.length; i < fd; i++) {
              FS.streams[i] = null; // Keep dense
            }
          }
        }
        // Close WebSocket first if we are about to replace the fd (i.e. dup2)
        if (FS.streams[fd] && FS.streams[fd].socket && FS.streams[fd].socket.close) {
          FS.streams[fd].socket.close();
        }
        FS.streams[fd] = stream;
        return fd;
      },removeFileHandle:function (fd) {
        FS.streams[fd] = null;
      },joinPath:function (parts, forceRelative) {
        var ret = parts[0];
        for (var i = 1; i < parts.length; i++) {
          if (ret[ret.length-1] != '/') ret += '/';
          ret += parts[i];
        }
        if (forceRelative && ret[0] == '/') ret = ret.substr(1);
        return ret;
      },absolutePath:function (relative, base) {
        if (typeof relative !== 'string') return null;
        if (base === undefined) base = FS.currentPath;
        if (relative && relative[0] == '/') base = '';
        var full = base + '/' + relative;
        var parts = full.split('/').reverse();
        var absolute = [''];
        while (parts.length) {
          var part = parts.pop();
          if (part == '' || part == '.') {
            // Nothing.
          } else if (part == '..') {
            if (absolute.length > 1) absolute.pop();
          } else {
            absolute.push(part);
          }
        }
        return absolute.length == 1 ? '/' : absolute.join('/');
      },analyzePath:function (path, dontResolveLastLink, linksVisited) {
        var ret = {
          isRoot: false,
          exists: false,
          error: 0,
          name: null,
          path: null,
          object: null,
          parentExists: false,
          parentPath: null,
          parentObject: null
        };
        path = FS.absolutePath(path);
        if (path == '/') {
          ret.isRoot = true;
          ret.exists = ret.parentExists = true;
          ret.name = '/';
          ret.path = ret.parentPath = '/';
          ret.object = ret.parentObject = FS.root;
        } else if (path !== null) {
          linksVisited = linksVisited || 0;
          path = path.slice(1).split('/');
          var current = FS.root;
          var traversed = [''];
          while (path.length) {
            if (path.length == 1 && current.isFolder) {
              ret.parentExists = true;
              ret.parentPath = traversed.length == 1 ? '/' : traversed.join('/');
              ret.parentObject = current;
              ret.name = path[0];
            }
            var target = path.shift();
            if (!current.isFolder) {
              ret.error = ERRNO_CODES.ENOTDIR;
              break;
            } else if (!current.read) {
              ret.error = ERRNO_CODES.EACCES;
              break;
            } else if (!current.contents.hasOwnProperty(target)) {
              ret.error = ERRNO_CODES.ENOENT;
              break;
            }
            current = current.contents[target];
            if (current.link && !(dontResolveLastLink && path.length == 0)) {
              if (linksVisited > 40) { // Usual Linux SYMLOOP_MAX.
                ret.error = ERRNO_CODES.ELOOP;
                break;
              }
              var link = FS.absolutePath(current.link, traversed.join('/'));
              ret = FS.analyzePath([link].concat(path).join('/'),
                                   dontResolveLastLink, linksVisited + 1);
              return ret;
            }
            traversed.push(target);
            if (path.length == 0) {
              ret.exists = true;
              ret.path = traversed.join('/');
              ret.object = current;
            }
          }
        }
        return ret;
      },findObject:function (path, dontResolveLastLink) {
        FS.ensureRoot();
        var ret = FS.analyzePath(path, dontResolveLastLink);
        if (ret.exists) {
          return ret.object;
        } else {
          ___setErrNo(ret.error);
          return null;
        }
      },createObject:function (parent, name, properties, canRead, canWrite) {
        if (!parent) parent = '/';
        if (typeof parent === 'string') parent = FS.findObject(parent);
        if (!parent) {
          ___setErrNo(ERRNO_CODES.EACCES);
          throw new Error('Parent path must exist.');
        }
        if (!parent.isFolder) {
          ___setErrNo(ERRNO_CODES.ENOTDIR);
          throw new Error('Parent must be a folder.');
        }
        if (!parent.write && !FS.ignorePermissions) {
          ___setErrNo(ERRNO_CODES.EACCES);
          throw new Error('Parent folder must be writeable.');
        }
        if (!name || name == '.' || name == '..') {
          ___setErrNo(ERRNO_CODES.ENOENT);
          throw new Error('Name must not be empty.');
        }
        if (parent.contents.hasOwnProperty(name)) {
          ___setErrNo(ERRNO_CODES.EEXIST);
          throw new Error("Can't overwrite object.");
        }
        parent.contents[name] = {
          read: canRead === undefined ? true : canRead,
          write: canWrite === undefined ? false : canWrite,
          timestamp: Date.now(),
          inodeNumber: FS.nextInode++
        };
        for (var key in properties) {
          if (properties.hasOwnProperty(key)) {
            parent.contents[name][key] = properties[key];
          }
        }
        return parent.contents[name];
      },createFolder:function (parent, name, canRead, canWrite) {
        var properties = {isFolder: true, isDevice: false, contents: {}};
        return FS.createObject(parent, name, properties, canRead, canWrite);
      },createPath:function (parent, path, canRead, canWrite) {
        var current = FS.findObject(parent);
        if (current === null) throw new Error('Invalid parent.');
        path = path.split('/').reverse();
        while (path.length) {
          var part = path.pop();
          if (!part) continue;
          if (!current.contents.hasOwnProperty(part)) {
            FS.createFolder(current, part, canRead, canWrite);
          }
          current = current.contents[part];
        }
        return current;
      },createFile:function (parent, name, properties, canRead, canWrite) {
        properties.isFolder = false;
        return FS.createObject(parent, name, properties, canRead, canWrite);
      },createDataFile:function (parent, name, data, canRead, canWrite) {
        if (typeof data === 'string') {
          var dataArray = new Array(data.length);
          for (var i = 0, len = data.length; i < len; ++i) dataArray[i] = data.charCodeAt(i);
          data = dataArray;
        }
        var properties = {
          isDevice: false,
          contents: data.subarray ? data.subarray(0) : data // as an optimization, create a new array wrapper (not buffer) here, to help JS engines understand this object
        };
        return FS.createFile(parent, name, properties, canRead, canWrite);
      },createLazyFile:function (parent, name, url, canRead, canWrite) {
        if (typeof XMLHttpRequest !== 'undefined') {
          if (!ENVIRONMENT_IS_WORKER) throw 'Cannot do synchronous binary XHRs outside webworkers in modern browsers. Use --embed-file or --preload-file in emcc';
          // Lazy chunked Uint8Array (implements get and length from Uint8Array). Actual getting is abstracted away for eventual reuse.
          var LazyUint8Array = function() {
            this.lengthKnown = false;
            this.chunks = []; // Loaded chunks. Index is the chunk number
          }
          LazyUint8Array.prototype.get = function(idx) {
            if (idx > this.length-1 || idx < 0) {
              return undefined;
            }
            var chunkOffset = idx % this.chunkSize;
            var chunkNum = Math.floor(idx / this.chunkSize);
            return this.getter(chunkNum)[chunkOffset];
          }
          LazyUint8Array.prototype.setDataGetter = function(getter) {
            this.getter = getter;
          }
          LazyUint8Array.prototype.cacheLength = function() {
              // Find length
              var xhr = new XMLHttpRequest();
              xhr.open('HEAD', url, false);
              xhr.send(null);
              if (!(xhr.status >= 200 && xhr.status < 300 || xhr.status === 304)) throw new Error("Couldn't load " + url + ". Status: " + xhr.status);
              var datalength = Number(xhr.getResponseHeader("Content-length"));
              var header;
              var hasByteServing = (header = xhr.getResponseHeader("Accept-Ranges")) && header === "bytes";
              var chunkSize = 1024*1024; // Chunk size in bytes
              if (!hasByteServing) chunkSize = datalength;
              // Function to get a range from the remote URL.
              var doXHR = (function(from, to) {
                if (from > to) throw new Error("invalid range (" + from + ", " + to + ") or no bytes requested!");
                if (to > datalength-1) throw new Error("only " + datalength + " bytes available! programmer error!");
                // TODO: Use mozResponseArrayBuffer, responseStream, etc. if available.
                var xhr = new XMLHttpRequest();
                xhr.open('GET', url, false);
                if (datalength !== chunkSize) xhr.setRequestHeader("Range", "bytes=" + from + "-" + to);
                // Some hints to the browser that we want binary data.
                if (typeof Uint8Array != 'undefined') xhr.responseType = 'arraybuffer';
                if (xhr.overrideMimeType) {
                  xhr.overrideMimeType('text/plain; charset=x-user-defined');
                }
                xhr.send(null);
                if (!(xhr.status >= 200 && xhr.status < 300 || xhr.status === 304)) throw new Error("Couldn't load " + url + ". Status: " + xhr.status);
                if (xhr.response !== undefined) {
                  return new Uint8Array(xhr.response || []);
                } else {
                  return intArrayFromString(xhr.responseText || '', true);
                }
              });
              var lazyArray = this;
              lazyArray.setDataGetter(function(chunkNum) {
                var start = chunkNum * chunkSize;
                var end = (chunkNum+1) * chunkSize - 1; // including this byte
                end = Math.min(end, datalength-1); // if datalength-1 is selected, this is the last block
                if (typeof(lazyArray.chunks[chunkNum]) === "undefined") {
                  lazyArray.chunks[chunkNum] = doXHR(start, end);
                }
                if (typeof(lazyArray.chunks[chunkNum]) === "undefined") throw new Error("doXHR failed!");
                return lazyArray.chunks[chunkNum];
              });
              this._length = datalength;
              this._chunkSize = chunkSize;
              this.lengthKnown = true;
          }
          var lazyArray = new LazyUint8Array();
          Object.defineProperty(lazyArray, "length", {
              get: function() {
                  if(!this.lengthKnown) {
                      this.cacheLength();
                  }
                  return this._length;
              }
          });
          Object.defineProperty(lazyArray, "chunkSize", {
              get: function() {
                  if(!this.lengthKnown) {
                      this.cacheLength();
                  }
                  return this._chunkSize;
              }
          });
          var properties = { isDevice: false, contents: lazyArray };
        } else {
          var properties = { isDevice: false, url: url };
        }
        return FS.createFile(parent, name, properties, canRead, canWrite);
      },createPreloadedFile:function (parent, name, url, canRead, canWrite, onload, onerror, dontCreateFile) {
        Browser.init();
        var fullname = FS.joinPath([parent, name], true);
        function processData(byteArray) {
          function finish(byteArray) {
            if (!dontCreateFile) {
              FS.createDataFile(parent, name, byteArray, canRead, canWrite);
            }
            if (onload) onload();
            removeRunDependency('cp ' + fullname);
          }
          var handled = false;
          Module['preloadPlugins'].forEach(function(plugin) {
            if (handled) return;
            if (plugin['canHandle'](fullname)) {
              plugin['handle'](byteArray, fullname, finish, function() {
                if (onerror) onerror();
                removeRunDependency('cp ' + fullname);
              });
              handled = true;
            }
          });
          if (!handled) finish(byteArray);
        }
        addRunDependency('cp ' + fullname);
        if (typeof url == 'string') {
          Browser.asyncLoad(url, function(byteArray) {
            processData(byteArray);
          }, onerror);
        } else {
          processData(url);
        }
      },createLink:function (parent, name, target, canRead, canWrite) {
        var properties = {isDevice: false, link: target};
        return FS.createFile(parent, name, properties, canRead, canWrite);
      },createDevice:function (parent, name, input, output) {
        if (!(input || output)) {
          throw new Error('A device must have at least one callback defined.');
        }
        var ops = {isDevice: true, input: input, output: output};
        return FS.createFile(parent, name, ops, Boolean(input), Boolean(output));
      },forceLoadFile:function (obj) {
        if (obj.isDevice || obj.isFolder || obj.link || obj.contents) return true;
        var success = true;
        if (typeof XMLHttpRequest !== 'undefined') {
          throw new Error("Lazy loading should have been performed (contents set) in createLazyFile, but it was not. Lazy loading only works in web workers. Use --embed-file or --preload-file in emcc on the main thread.");
        } else if (Module['read']) {
          // Command-line.
          try {
            // WARNING: Can't read binary files in V8's d8 or tracemonkey's js, as
            //          read() will try to parse UTF8.
            obj.contents = intArrayFromString(Module['read'](obj.url), true);
          } catch (e) {
            success = false;
          }
        } else {
          throw new Error('Cannot load without read() or XMLHttpRequest.');
        }
        if (!success) ___setErrNo(ERRNO_CODES.EIO);
        return success;
      },ensureRoot:function () {
        if (FS.root) return;
        // The main file system tree. All the contents are inside this.
        FS.root = {
          read: true,
          write: true,
          isFolder: true,
          isDevice: false,
          timestamp: Date.now(),
          inodeNumber: 1,
          contents: {}
        };
      },init:function (input, output, error) {
        // Make sure we initialize only once.
        assert(!FS.init.initialized, 'FS.init was previously called. If you want to initialize later with custom parameters, remove any earlier calls (note that one is automatically added to the generated code)');
        FS.init.initialized = true;
        FS.ensureRoot();
        // Allow Module.stdin etc. to provide defaults, if none explicitly passed to us here
        input = input || Module['stdin'];
        output = output || Module['stdout'];
        error = error || Module['stderr'];
        // Default handlers.
        var stdinOverridden = true, stdoutOverridden = true, stderrOverridden = true;
        if (!input) {
          stdinOverridden = false;
          input = function() {
            if (!input.cache || !input.cache.length) {
              var result;
              if (typeof window != 'undefined' &&
                  typeof window.prompt == 'function') {
                // Browser.
                result = window.prompt('Input: ');
                if (result === null) result = String.fromCharCode(0); // cancel ==> EOF
              } else if (typeof readline == 'function') {
                // Command line.
                result = readline();
              }
              if (!result) result = '';
              input.cache = intArrayFromString(result + '\n', true);
            }
            return input.cache.shift();
          };
        }
        var utf8 = new Runtime.UTF8Processor();
        function simpleOutput(val) {
          if (val === null || val === 10) {
            output.printer(output.buffer.join(''));
            output.buffer = [];
          } else {
            output.buffer.push(utf8.processCChar(val));
          }
        }
        if (!output) {
          stdoutOverridden = false;
          output = simpleOutput;
        }
        if (!output.printer) output.printer = Module['print'];
        if (!output.buffer) output.buffer = [];
        if (!error) {
          stderrOverridden = false;
          error = simpleOutput;
        }
        if (!error.printer) error.printer = Module['print'];
        if (!error.buffer) error.buffer = [];
        // Create the temporary folder, if not already created
        try {
          FS.createFolder('/', 'tmp', true, true);
        } catch(e) {}
        // Create the I/O devices.
        var devFolder = FS.createFolder('/', 'dev', true, true);
        var stdin = FS.createDevice(devFolder, 'stdin', input);
        var stdout = FS.createDevice(devFolder, 'stdout', null, output);
        var stderr = FS.createDevice(devFolder, 'stderr', null, error);
        FS.createDevice(devFolder, 'tty', input, output);
        FS.createDevice(devFolder, 'null', function(){}, function(){});
        // Create default streams.
        FS.streams[1] = {
          path: '/dev/stdin',
          object: stdin,
          position: 0,
          isRead: true,
          isWrite: false,
          isAppend: false,
          isTerminal: !stdinOverridden,
          error: false,
          eof: false,
          ungotten: []
        };
        FS.streams[2] = {
          path: '/dev/stdout',
          object: stdout,
          position: 0,
          isRead: false,
          isWrite: true,
          isAppend: false,
          isTerminal: !stdoutOverridden,
          error: false,
          eof: false,
          ungotten: []
        };
        FS.streams[3] = {
          path: '/dev/stderr',
          object: stderr,
          position: 0,
          isRead: false,
          isWrite: true,
          isAppend: false,
          isTerminal: !stderrOverridden,
          error: false,
          eof: false,
          ungotten: []
        };
        // TODO: put these low in memory like we used to assert on: assert(Math.max(_stdin, _stdout, _stderr) < 15000); // make sure these are low, we flatten arrays with these
        HEAP32[((_stdin)>>2)]=1;
        HEAP32[((_stdout)>>2)]=2;
        HEAP32[((_stderr)>>2)]=3;
        // Other system paths
        FS.createPath('/', 'dev/shm/tmp', true, true); // temp files
        // Newlib initialization
        for (var i = FS.streams.length; i < Math.max(_stdin, _stdout, _stderr) + 4; i++) {
          FS.streams[i] = null; // Make sure to keep FS.streams dense
        }
        FS.streams[_stdin] = FS.streams[1];
        FS.streams[_stdout] = FS.streams[2];
        FS.streams[_stderr] = FS.streams[3];
        allocate([ allocate(
          [0, 0, 0, 0, _stdin, 0, 0, 0, _stdout, 0, 0, 0, _stderr, 0, 0, 0],
          'void*', ALLOC_NORMAL) ], 'void*', ALLOC_NONE, __impure_ptr);
      },quit:function () {
        if (!FS.init.initialized) return;
        // Flush any partially-printed lines in stdout and stderr. Careful, they may have been closed
        if (FS.streams[2] && FS.streams[2].object.output.buffer.length > 0) FS.streams[2].object.output(10);
        if (FS.streams[3] && FS.streams[3].object.output.buffer.length > 0) FS.streams[3].object.output(10);
      },standardizePath:function (path) {
        if (path.substr(0, 2) == './') path = path.substr(2);
        return path;
      },deleteFile:function (path) {
        path = FS.analyzePath(path);
        if (!path.parentExists || !path.exists) {
          throw 'Invalid path ' + path;
        }
        delete path.parentObject.contents[path.name];
      }};
  var ___dirent_struct_layout={__size__:1040,d_ino:0,d_name:4,d_off:1028,d_reclen:1032,d_type:1036};function _open(path, oflag, varargs) {
      // int open(const char *path, int oflag, ...);
      // http://pubs.opengroup.org/onlinepubs/009695399/functions/open.html
      // NOTE: This implementation tries to mimic glibc rather than strictly
      // following the POSIX standard.
      var mode = HEAP32[((varargs)>>2)];
      // Simplify flags.
      var accessMode = oflag & 3;
      var isWrite = accessMode != 0;
      var isRead = accessMode != 1;
      var isCreate = Boolean(oflag & 512);
      var isExistCheck = Boolean(oflag & 2048);
      var isTruncate = Boolean(oflag & 1024);
      var isAppend = Boolean(oflag & 8);
      // Verify path.
      var origPath = path;
      path = FS.analyzePath(Pointer_stringify(path));
      if (!path.parentExists) {
        ___setErrNo(path.error);
        return -1;
      }
      var target = path.object || null;
      var finalPath;
      // Verify the file exists, create if needed and allowed.
      if (target) {
        if (isCreate && isExistCheck) {
          ___setErrNo(ERRNO_CODES.EEXIST);
          return -1;
        }
        if ((isWrite || isCreate || isTruncate) && target.isFolder) {
          ___setErrNo(ERRNO_CODES.EISDIR);
          return -1;
        }
        if (isRead && !target.read || isWrite && !target.write) {
          ___setErrNo(ERRNO_CODES.EACCES);
          return -1;
        }
        if (isTruncate && !target.isDevice) {
          target.contents = [];
        } else {
          if (!FS.forceLoadFile(target)) {
            ___setErrNo(ERRNO_CODES.EIO);
            return -1;
          }
        }
        finalPath = path.path;
      } else {
        if (!isCreate) {
          ___setErrNo(ERRNO_CODES.ENOENT);
          return -1;
        }
        if (!path.parentObject.write) {
          ___setErrNo(ERRNO_CODES.EACCES);
          return -1;
        }
        target = FS.createDataFile(path.parentObject, path.name, [],
                                   mode & 0x100, mode & 0x80);  // S_IRUSR, S_IWUSR.
        finalPath = path.parentPath + '/' + path.name;
      }
      // Actually create an open stream.
      var id;
      if (target.isFolder) {
        var entryBuffer = 0;
        if (___dirent_struct_layout) {
          entryBuffer = _malloc(___dirent_struct_layout.__size__);
        }
        var contents = [];
        for (var key in target.contents) contents.push(key);
        id = FS.createFileHandle({
          path: finalPath,
          object: target,
          // An index into contents. Special values: -2 is ".", -1 is "..".
          position: -2,
          isRead: true,
          isWrite: false,
          isAppend: false,
          error: false,
          eof: false,
          ungotten: [],
          // Folder-specific properties:
          // Remember the contents at the time of opening in an array, so we can
          // seek between them relying on a single order.
          contents: contents,
          // Each stream has its own area for readdir() returns.
          currentEntry: entryBuffer
        });
      } else {
        id = FS.createFileHandle({
          path: finalPath,
          object: target,
          position: 0,
          isRead: isRead,
          isWrite: isWrite,
          isAppend: isAppend,
          error: false,
          eof: false,
          ungotten: []
        });
      }
      return id;
    }function _fopen(filename, mode) {
      // FILE *fopen(const char *restrict filename, const char *restrict mode);
      // http://pubs.opengroup.org/onlinepubs/000095399/functions/fopen.html
      var flags;
      mode = Pointer_stringify(mode);
      if (mode[0] == 'r') {
        if (mode.indexOf('+') != -1) {
          flags = 2;
        } else {
          flags = 0;
        }
      } else if (mode[0] == 'w') {
        if (mode.indexOf('+') != -1) {
          flags = 2;
        } else {
          flags = 1;
        }
        flags |= 512;
        flags |= 1024;
      } else if (mode[0] == 'a') {
        if (mode.indexOf('+') != -1) {
          flags = 2;
        } else {
          flags = 1;
        }
        flags |= 512;
        flags |= 8;
      } else {
        ___setErrNo(ERRNO_CODES.EINVAL);
        return 0;
      }
      var ret = _open(filename, flags, allocate([0x1FF, 0, 0, 0], 'i32', ALLOC_STACK));  // All creation permissions.
      return (ret == -1) ? 0 : ret;
    }
  function _close(fildes) {
      // int close(int fildes);
      // http://pubs.opengroup.org/onlinepubs/000095399/functions/close.html
      if (FS.streams[fildes]) {
        if (FS.streams[fildes].currentEntry) {
          _free(FS.streams[fildes].currentEntry);
        }
        FS.streams[fildes] = null;
        return 0;
      } else {
        ___setErrNo(ERRNO_CODES.EBADF);
        return -1;
      }
    }
  function _fsync(fildes) {
      // int fsync(int fildes);
      // http://pubs.opengroup.org/onlinepubs/000095399/functions/fsync.html
      if (FS.streams[fildes]) {
        // We write directly to the file system, so there's nothing to do here.
        return 0;
      } else {
        ___setErrNo(ERRNO_CODES.EBADF);
        return -1;
      }
    }function _fclose(stream) {
      // int fclose(FILE *stream);
      // http://pubs.opengroup.org/onlinepubs/000095399/functions/fclose.html
      _fsync(stream);
      return _close(stream);
    }
;
;
  function __exit(status) {
      // void _exit(int status);
      // http://pubs.opengroup.org/onlinepubs/000095399/functions/exit.html
      function ExitStatus() {
        this.name = "ExitStatus";
        this.message = "Program terminated with exit(" + status + ")";
        this.status = status;
        Module.print('Exit Status: ' + status);
      };
      ExitStatus.prototype = new Error();
      ExitStatus.prototype.constructor = ExitStatus;
      exitRuntime();
      ABORT = true;
      throw new ExitStatus();
    }function _exit(status) {
      __exit(status);
    }function __ZSt9terminatev() {
      _exit(-1234);
    }
;
  function _memset(ptr, value, num) {
      ptr = ptr|0; value = value|0; num = num|0;
      var stop = 0, value4 = 0, stop4 = 0, unaligned = 0;
      stop = (ptr + num)|0;
      if ((num|0) >= 20) {
        // This is unaligned, but quite large, so work hard to get to aligned settings
        value = value & 0xff;
        unaligned = ptr & 3;
        value4 = value | (value << 8) | (value << 16) | (value << 24);
        stop4 = stop & ~3;
        if (unaligned) {
          unaligned = (ptr + 4 - unaligned)|0;
          while ((ptr|0) < (unaligned|0)) { // no need to check for stop, since we have large num
            HEAP8[(ptr)]=value;
            ptr = (ptr+1)|0;
          }
        }
        while ((ptr|0) < (stop4|0)) {
          HEAP32[((ptr)>>2)]=value4;
          ptr = (ptr+4)|0;
        }
      }
      while ((ptr|0) < (stop|0)) {
        HEAP8[(ptr)]=value;
        ptr = (ptr+1)|0;
      }
    }var _llvm_memset_p0i8_i32=_memset;
;
  function ___assert_func(filename, line, func, condition) {
      throw 'Assertion failed: ' + (condition ? Pointer_stringify(condition) : 'unknown condition') + ', at: ' + [filename ? Pointer_stringify(filename) : 'unknown filename', line, func ? Pointer_stringify(func) : 'unknown function'] + ' at ' + new Error().stack;
    }
  function _strlen(ptr) {
      ptr = ptr|0;
      var curr = 0;
      curr = ptr;
      while (HEAP8[(curr)]) {
        curr = (curr + 1)|0;
      }
      return (curr - ptr)|0;
    }
  var _llvm_memcpy_p0i8_p0i8_i64=_memcpy;
  function _fflush(stream) {
      // int fflush(FILE *stream);
      // http://pubs.opengroup.org/onlinepubs/000095399/functions/fflush.html
      var flush = function(filedes) {
        // Right now we write all data directly, except for output devices.
        if (FS.streams[filedes] && FS.streams[filedes].object.output) {
          if (!FS.streams[filedes].isTerminal) { // don't flush terminals, it would cause a \n to also appear
            FS.streams[filedes].object.output(null);
          }
        }
      };
      try {
        if (stream === 0) {
          for (var i = 0; i < FS.streams.length; i++) if (FS.streams[i]) flush(i);
        } else {
          flush(stream);
        }
        return 0;
      } catch (e) {
        ___setErrNo(ERRNO_CODES.EIO);
        return -1;
      }
    }
  function _lseek(fildes, offset, whence) {
      // off_t lseek(int fildes, off_t offset, int whence);
      // http://pubs.opengroup.org/onlinepubs/000095399/functions/lseek.html
      if (FS.streams[fildes] && !FS.streams[fildes].object.isDevice) {
        var stream = FS.streams[fildes];
        var position = offset;
        if (whence === 1) {  // SEEK_CUR.
          position += stream.position;
        } else if (whence === 2) {  // SEEK_END.
          position += stream.object.contents.length;
        }
        if (position < 0) {
          ___setErrNo(ERRNO_CODES.EINVAL);
          return -1;
        } else {
          stream.ungotten = [];
          stream.position = position;
          return position;
        }
      } else {
        ___setErrNo(ERRNO_CODES.EBADF);
        return -1;
      }
    }function _fseek(stream, offset, whence) {
      // int fseek(FILE *stream, long offset, int whence);
      // http://pubs.opengroup.org/onlinepubs/000095399/functions/fseek.html
      var ret = _lseek(stream, offset, whence);
      if (ret == -1) {
        return -1;
      } else {
        FS.streams[stream].eof = false;
        return 0;
      }
    }
  function _recv(fd, buf, len, flags) {
      var info = FS.streams[fd];
      if (!info) return -1;
      if (!info.hasData()) {
        ___setErrNo(ERRNO_CODES.EAGAIN); // no data, and all sockets are nonblocking, so this is the right behavior
        return -1;
      }
      var buffer = info.inQueue.shift();
      if (len < buffer.length) {
        if (info.stream) {
          // This is tcp (reliable), so if not all was read, keep it
          info.inQueue.unshift(buffer.subarray(len));
        }
        buffer = buffer.subarray(0, len);
      }
      HEAPU8.set(buffer, buf);
      return buffer.length;
    }
  function _pread(fildes, buf, nbyte, offset) {
      // ssize_t pread(int fildes, void *buf, size_t nbyte, off_t offset);
      // http://pubs.opengroup.org/onlinepubs/000095399/functions/read.html
      var stream = FS.streams[fildes];
      if (!stream || stream.object.isDevice) {
        ___setErrNo(ERRNO_CODES.EBADF);
        return -1;
      } else if (!stream.isRead) {
        ___setErrNo(ERRNO_CODES.EACCES);
        return -1;
      } else if (stream.object.isFolder) {
        ___setErrNo(ERRNO_CODES.EISDIR);
        return -1;
      } else if (nbyte < 0 || offset < 0) {
        ___setErrNo(ERRNO_CODES.EINVAL);
        return -1;
      } else {
        var bytesRead = 0;
        while (stream.ungotten.length && nbyte > 0) {
          HEAP8[((buf++)|0)]=stream.ungotten.pop()
          nbyte--;
          bytesRead++;
        }
        var contents = stream.object.contents;
        var size = Math.min(contents.length - offset, nbyte);
        if (contents.subarray) { // typed array
          HEAPU8.set(contents.subarray(offset, offset+size), buf);
        } else
        if (contents.slice) { // normal array
          for (var i = 0; i < size; i++) {
            HEAP8[(((buf)+(i))|0)]=contents[offset + i]
          }
        } else {
          for (var i = 0; i < size; i++) { // LazyUint8Array from sync binary XHR
            HEAP8[(((buf)+(i))|0)]=contents.get(offset + i)
          }
        }
        bytesRead += size;
        return bytesRead;
      }
    }function _read(fildes, buf, nbyte) {
      // ssize_t read(int fildes, void *buf, size_t nbyte);
      // http://pubs.opengroup.org/onlinepubs/000095399/functions/read.html
      var stream = FS.streams[fildes];
      if (stream && ('socket' in stream)) {
        return _recv(fildes, buf, nbyte, 0);
      } else if (!stream) {
        ___setErrNo(ERRNO_CODES.EBADF);
        return -1;
      } else if (!stream.isRead) {
        ___setErrNo(ERRNO_CODES.EACCES);
        return -1;
      } else if (nbyte < 0) {
        ___setErrNo(ERRNO_CODES.EINVAL);
        return -1;
      } else {
        var bytesRead;
        if (stream.object.isDevice) {
          if (stream.object.input) {
            bytesRead = 0;
            while (stream.ungotten.length && nbyte > 0) {
              HEAP8[((buf++)|0)]=stream.ungotten.pop()
              nbyte--;
              bytesRead++;
            }
            for (var i = 0; i < nbyte; i++) {
              try {
                var result = stream.object.input();
              } catch (e) {
                ___setErrNo(ERRNO_CODES.EIO);
                return -1;
              }
              if (result === undefined && bytesRead === 0) {
                ___setErrNo(ERRNO_CODES.EAGAIN);
                return -1;
              }
              if (result === null || result === undefined) break;
              bytesRead++;
              HEAP8[(((buf)+(i))|0)]=result
            }
            return bytesRead;
          } else {
            ___setErrNo(ERRNO_CODES.ENXIO);
            return -1;
          }
        } else {
          var ungotSize = stream.ungotten.length;
          bytesRead = _pread(fildes, buf, nbyte, stream.position);
          if (bytesRead != -1) {
            stream.position += (stream.ungotten.length - ungotSize) + bytesRead;
          }
          return bytesRead;
        }
      }
    }function _fread(ptr, size, nitems, stream) {
      // size_t fread(void *restrict ptr, size_t size, size_t nitems, FILE *restrict stream);
      // http://pubs.opengroup.org/onlinepubs/000095399/functions/fread.html
      var bytesToRead = nitems * size;
      if (bytesToRead == 0) return 0;
      var bytesRead = _read(stream, ptr, bytesToRead);
      var streamObj = FS.streams[stream];
      if (bytesRead == -1) {
        if (streamObj) streamObj.error = true;
        return 0;
      } else {
        if (bytesRead < bytesToRead) streamObj.eof = true;
        return Math.floor(bytesRead / size);
      }
    }
  var _llvm_va_start=undefined;
  function _send(fd, buf, len, flags) {
      var info = FS.streams[fd];
      if (!info) return -1;
      info.sender(HEAPU8.subarray(buf, buf+len));
      return len;
    }
  function _pwrite(fildes, buf, nbyte, offset) {
      // ssize_t pwrite(int fildes, const void *buf, size_t nbyte, off_t offset);
      // http://pubs.opengroup.org/onlinepubs/000095399/functions/write.html
      var stream = FS.streams[fildes];
      if (!stream || stream.object.isDevice) {
        ___setErrNo(ERRNO_CODES.EBADF);
        return -1;
      } else if (!stream.isWrite) {
        ___setErrNo(ERRNO_CODES.EACCES);
        return -1;
      } else if (stream.object.isFolder) {
        ___setErrNo(ERRNO_CODES.EISDIR);
        return -1;
      } else if (nbyte < 0 || offset < 0) {
        ___setErrNo(ERRNO_CODES.EINVAL);
        return -1;
      } else {
        var contents = stream.object.contents;
        while (contents.length < offset) contents.push(0);
        for (var i = 0; i < nbyte; i++) {
          contents[offset + i] = HEAPU8[(((buf)+(i))|0)];
        }
        stream.object.timestamp = Date.now();
        return i;
      }
    }function _write(fildes, buf, nbyte) {
      // ssize_t write(int fildes, const void *buf, size_t nbyte);
      // http://pubs.opengroup.org/onlinepubs/000095399/functions/write.html
      var stream = FS.streams[fildes];
      if (stream && ('socket' in stream)) {
          return _send(fildes, buf, nbyte, 0);
      } else if (!stream) {
        ___setErrNo(ERRNO_CODES.EBADF);
        return -1;
      } else if (!stream.isWrite) {
        ___setErrNo(ERRNO_CODES.EACCES);
        return -1;
      } else if (nbyte < 0) {
        ___setErrNo(ERRNO_CODES.EINVAL);
        return -1;
      } else {
        if (stream.object.isDevice) {
          if (stream.object.output) {
            for (var i = 0; i < nbyte; i++) {
              try {
                stream.object.output(HEAP8[(((buf)+(i))|0)]);
              } catch (e) {
                ___setErrNo(ERRNO_CODES.EIO);
                return -1;
              }
            }
            stream.object.timestamp = Date.now();
            return i;
          } else {
            ___setErrNo(ERRNO_CODES.ENXIO);
            return -1;
          }
        } else {
          var bytesWritten = _pwrite(fildes, buf, nbyte, stream.position);
          if (bytesWritten != -1) stream.position += bytesWritten;
          return bytesWritten;
        }
      }
    }function _fwrite(ptr, size, nitems, stream) {
      // size_t fwrite(const void *restrict ptr, size_t size, size_t nitems, FILE *restrict stream);
      // http://pubs.opengroup.org/onlinepubs/000095399/functions/fwrite.html
      var bytesToWrite = nitems * size;
      if (bytesToWrite == 0) return 0;
      var bytesWritten = _write(stream, ptr, bytesToWrite);
      if (bytesWritten == -1) {
        if (FS.streams[stream]) FS.streams[stream].error = true;
        return 0;
      } else {
        return Math.floor(bytesWritten / size);
      }
    }
  function __reallyNegative(x) {
      return x < 0 || (x === 0 && (1/x) === -Infinity);
    }function __formatString(format, varargs) {
      var textIndex = format;
      var argIndex = 0;
      function getNextArg(type) {
        // NOTE: Explicitly ignoring type safety. Otherwise this fails:
        //       int x = 4; printf("%c\n", (char)x);
        var ret;
        if (type === 'double') {
          ret = HEAPF64[(((varargs)+(argIndex))>>3)];
        } else if (type == 'i64') {
          ret = [HEAP32[(((varargs)+(argIndex))>>2)],
                 HEAP32[(((varargs)+(argIndex+8))>>2)]];
          argIndex += 8; // each 32-bit chunk is in a 64-bit block
        } else {
          type = 'i32'; // varargs are always i32, i64, or double
          ret = HEAP32[(((varargs)+(argIndex))>>2)];
        }
        argIndex += Math.max(Runtime.getNativeFieldSize(type), Runtime.getAlignSize(type, null, true));
        return ret;
      }
      var ret = [];
      var curr, next, currArg;
      while(1) {
        var startTextIndex = textIndex;
        curr = HEAP8[(textIndex)];
        if (curr === 0) break;
        next = HEAP8[((textIndex+1)|0)];
        if (curr == 37) {
          // Handle flags.
          var flagAlwaysSigned = false;
          var flagLeftAlign = false;
          var flagAlternative = false;
          var flagZeroPad = false;
          flagsLoop: while (1) {
            switch (next) {
              case 43:
                flagAlwaysSigned = true;
                break;
              case 45:
                flagLeftAlign = true;
                break;
              case 35:
                flagAlternative = true;
                break;
              case 48:
                if (flagZeroPad) {
                  break flagsLoop;
                } else {
                  flagZeroPad = true;
                  break;
                }
              default:
                break flagsLoop;
            }
            textIndex++;
            next = HEAP8[((textIndex+1)|0)];
          }
          // Handle width.
          var width = 0;
          if (next == 42) {
            width = getNextArg('i32');
            textIndex++;
            next = HEAP8[((textIndex+1)|0)];
          } else {
            while (next >= 48 && next <= 57) {
              width = width * 10 + (next - 48);
              textIndex++;
              next = HEAP8[((textIndex+1)|0)];
            }
          }
          // Handle precision.
          var precisionSet = false;
          if (next == 46) {
            var precision = 0;
            precisionSet = true;
            textIndex++;
            next = HEAP8[((textIndex+1)|0)];
            if (next == 42) {
              precision = getNextArg('i32');
              textIndex++;
            } else {
              while(1) {
                var precisionChr = HEAP8[((textIndex+1)|0)];
                if (precisionChr < 48 ||
                    precisionChr > 57) break;
                precision = precision * 10 + (precisionChr - 48);
                textIndex++;
              }
            }
            next = HEAP8[((textIndex+1)|0)];
          } else {
            var precision = 6; // Standard default.
          }
          // Handle integer sizes. WARNING: These assume a 32-bit architecture!
          var argSize;
          switch (String.fromCharCode(next)) {
            case 'h':
              var nextNext = HEAP8[((textIndex+2)|0)];
              if (nextNext == 104) {
                textIndex++;
                argSize = 1; // char (actually i32 in varargs)
              } else {
                argSize = 2; // short (actually i32 in varargs)
              }
              break;
            case 'l':
              var nextNext = HEAP8[((textIndex+2)|0)];
              if (nextNext == 108) {
                textIndex++;
                argSize = 8; // long long
              } else {
                argSize = 4; // long
              }
              break;
            case 'L': // long long
            case 'q': // int64_t
            case 'j': // intmax_t
              argSize = 8;
              break;
            case 'z': // size_t
            case 't': // ptrdiff_t
            case 'I': // signed ptrdiff_t or unsigned size_t
              argSize = 4;
              break;
            default:
              argSize = null;
          }
          if (argSize) textIndex++;
          next = HEAP8[((textIndex+1)|0)];
          // Handle type specifier.
          switch (String.fromCharCode(next)) {
            case 'd': case 'i': case 'u': case 'o': case 'x': case 'X': case 'p': {
              // Integer.
              var signed = next == 100 || next == 105;
              argSize = argSize || 4;
              var currArg = getNextArg('i' + (argSize * 8));
              var argText;
              // Flatten i64-1 [low, high] into a (slightly rounded) double
              if (argSize == 8) {
                currArg = Runtime.makeBigInt(currArg[0], currArg[1], next == 117);
              }
              // Truncate to requested size.
              if (argSize <= 4) {
                var limit = Math.pow(256, argSize) - 1;
                currArg = (signed ? reSign : unSign)(currArg & limit, argSize * 8);
              }
              // Format the number.
              var currAbsArg = Math.abs(currArg);
              var prefix = '';
              if (next == 100 || next == 105) {
                argText = reSign(currArg, 8 * argSize, 1).toString(10);
              } else if (next == 117) {
                argText = unSign(currArg, 8 * argSize, 1).toString(10);
                currArg = Math.abs(currArg);
              } else if (next == 111) {
                argText = (flagAlternative ? '0' : '') + currAbsArg.toString(8);
              } else if (next == 120 || next == 88) {
                prefix = (flagAlternative && currArg != 0) ? '0x' : '';
                if (currArg < 0) {
                  // Represent negative numbers in hex as 2's complement.
                  currArg = -currArg;
                  argText = (currAbsArg - 1).toString(16);
                  var buffer = [];
                  for (var i = 0; i < argText.length; i++) {
                    buffer.push((0xF - parseInt(argText[i], 16)).toString(16));
                  }
                  argText = buffer.join('');
                  while (argText.length < argSize * 2) argText = 'f' + argText;
                } else {
                  argText = currAbsArg.toString(16);
                }
                if (next == 88) {
                  prefix = prefix.toUpperCase();
                  argText = argText.toUpperCase();
                }
              } else if (next == 112) {
                if (currAbsArg === 0) {
                  argText = '(nil)';
                } else {
                  prefix = '0x';
                  argText = currAbsArg.toString(16);
                }
              }
              if (precisionSet) {
                while (argText.length < precision) {
                  argText = '0' + argText;
                }
              }
              // Add sign if needed
              if (flagAlwaysSigned) {
                if (currArg < 0) {
                  prefix = '-' + prefix;
                } else {
                  prefix = '+' + prefix;
                }
              }
              // Add padding.
              while (prefix.length + argText.length < width) {
                if (flagLeftAlign) {
                  argText += ' ';
                } else {
                  if (flagZeroPad) {
                    argText = '0' + argText;
                  } else {
                    prefix = ' ' + prefix;
                  }
                }
              }
              // Insert the result into the buffer.
              argText = prefix + argText;
              argText.split('').forEach(function(chr) {
                ret.push(chr.charCodeAt(0));
              });
              break;
            }
            case 'f': case 'F': case 'e': case 'E': case 'g': case 'G': {
              // Float.
              var currArg = getNextArg('double');
              var argText;
              if (isNaN(currArg)) {
                argText = 'nan';
                flagZeroPad = false;
              } else if (!isFinite(currArg)) {
                argText = (currArg < 0 ? '-' : '') + 'inf';
                flagZeroPad = false;
              } else {
                var isGeneral = false;
                var effectivePrecision = Math.min(precision, 20);
                // Convert g/G to f/F or e/E, as per:
                // http://pubs.opengroup.org/onlinepubs/9699919799/functions/printf.html
                if (next == 103 || next == 71) {
                  isGeneral = true;
                  precision = precision || 1;
                  var exponent = parseInt(currArg.toExponential(effectivePrecision).split('e')[1], 10);
                  if (precision > exponent && exponent >= -4) {
                    next = ((next == 103) ? 'f' : 'F').charCodeAt(0);
                    precision -= exponent + 1;
                  } else {
                    next = ((next == 103) ? 'e' : 'E').charCodeAt(0);
                    precision--;
                  }
                  effectivePrecision = Math.min(precision, 20);
                }
                if (next == 101 || next == 69) {
                  argText = currArg.toExponential(effectivePrecision);
                  // Make sure the exponent has at least 2 digits.
                  if (/[eE][-+]\d$/.test(argText)) {
                    argText = argText.slice(0, -1) + '0' + argText.slice(-1);
                  }
                } else if (next == 102 || next == 70) {
                  argText = currArg.toFixed(effectivePrecision);
                  if (currArg === 0 && __reallyNegative(currArg)) {
                    argText = '-' + argText;
                  }
                }
                var parts = argText.split('e');
                if (isGeneral && !flagAlternative) {
                  // Discard trailing zeros and periods.
                  while (parts[0].length > 1 && parts[0].indexOf('.') != -1 &&
                         (parts[0].slice(-1) == '0' || parts[0].slice(-1) == '.')) {
                    parts[0] = parts[0].slice(0, -1);
                  }
                } else {
                  // Make sure we have a period in alternative mode.
                  if (flagAlternative && argText.indexOf('.') == -1) parts[0] += '.';
                  // Zero pad until required precision.
                  while (precision > effectivePrecision++) parts[0] += '0';
                }
                argText = parts[0] + (parts.length > 1 ? 'e' + parts[1] : '');
                // Capitalize 'E' if needed.
                if (next == 69) argText = argText.toUpperCase();
                // Add sign.
                if (flagAlwaysSigned && currArg >= 0) {
                  argText = '+' + argText;
                }
              }
              // Add padding.
              while (argText.length < width) {
                if (flagLeftAlign) {
                  argText += ' ';
                } else {
                  if (flagZeroPad && (argText[0] == '-' || argText[0] == '+')) {
                    argText = argText[0] + '0' + argText.slice(1);
                  } else {
                    argText = (flagZeroPad ? '0' : ' ') + argText;
                  }
                }
              }
              // Adjust case.
              if (next < 97) argText = argText.toUpperCase();
              // Insert the result into the buffer.
              argText.split('').forEach(function(chr) {
                ret.push(chr.charCodeAt(0));
              });
              break;
            }
            case 's': {
              // String.
              var arg = getNextArg('i8*');
              var argLength = arg ? _strlen(arg) : '(null)'.length;
              if (precisionSet) argLength = Math.min(argLength, precision);
              if (!flagLeftAlign) {
                while (argLength < width--) {
                  ret.push(32);
                }
              }
              if (arg) {
                for (var i = 0; i < argLength; i++) {
                  ret.push(HEAPU8[((arg++)|0)]);
                }
              } else {
                ret = ret.concat(intArrayFromString('(null)'.substr(0, argLength), true));
              }
              if (flagLeftAlign) {
                while (argLength < width--) {
                  ret.push(32);
                }
              }
              break;
            }
            case 'c': {
              // Character.
              if (flagLeftAlign) ret.push(getNextArg('i8'));
              while (--width > 0) {
                ret.push(32);
              }
              if (!flagLeftAlign) ret.push(getNextArg('i8'));
              break;
            }
            case 'n': {
              // Write the length written so far to the next parameter.
              var ptr = getNextArg('i32*');
              HEAP32[((ptr)>>2)]=ret.length
              break;
            }
            case '%': {
              // Literal percent sign.
              ret.push(curr);
              break;
            }
            default: {
              // Unknown specifiers remain untouched.
              for (var i = startTextIndex; i < textIndex + 2; i++) {
                ret.push(HEAP8[(i)]);
              }
            }
          }
          textIndex += 2;
          // TODO: Support a/A (hex float) and m (last error) specifiers.
          // TODO: Support %1${specifier} for arg selection.
        } else {
          ret.push(curr);
          textIndex += 1;
        }
      }
      return ret;
    }function _fprintf(stream, format, varargs) {
      // int fprintf(FILE *restrict stream, const char *restrict format, ...);
      // http://pubs.opengroup.org/onlinepubs/000095399/functions/printf.html
      var result = __formatString(format, varargs);
      var stack = Runtime.stackSave();
      var ret = _fwrite(allocate(result, 'i8', ALLOC_STACK), 1, result.length, stream);
      Runtime.stackRestore(stack);
      return ret;
    }
  function _vfprintf(s, f, va_arg) {
      return _fprintf(s, f, HEAP32[((va_arg)>>2)]);
    }
  function _llvm_va_end() {}
  function _fputs(s, stream) {
      // int fputs(const char *restrict s, FILE *restrict stream);
      // http://pubs.opengroup.org/onlinepubs/000095399/functions/fputs.html
      return _write(stream, s, _strlen(s));
    }
  function _snprintf(s, n, format, varargs) {
      // int snprintf(char *restrict s, size_t n, const char *restrict format, ...);
      // http://pubs.opengroup.org/onlinepubs/000095399/functions/printf.html
      var result = __formatString(format, varargs);
      var limit = (n === undefined) ? result.length
                                    : Math.min(result.length, Math.max(n - 1, 0));
      if (s < 0) {
        s = -s;
        var buf = _malloc(limit+1);
        HEAP32[((s)>>2)]=buf;
        s = buf;
      }
      for (var i = 0; i < limit; i++) {
        HEAP8[(((s)+(i))|0)]=result[i];
      }
      if (limit < n || (n === undefined)) HEAP8[(((s)+(i))|0)]=0;
      return result.length;
    }function _vsnprintf(s, n, format, va_arg) {
      return _snprintf(s, n, format, HEAP32[((va_arg)>>2)]);
    }
  function _strcpy(pdest, psrc) {
      pdest = pdest|0; psrc = psrc|0;
      var i = 0;
      do {
        HEAP8[(((pdest+i)|0)|0)]=HEAP8[(((psrc+i)|0)|0)];
        i = (i+1)|0;
      } while (HEAP8[(((psrc)+(i-1))|0)]);
      return pdest|0;
    }
  function _strcat(pdest, psrc) {
      pdest = pdest|0; psrc = psrc|0;
      var i = 0;
      var pdestEnd = 0;
      pdestEnd = (pdest + (_strlen(pdest)|0))|0;
      do {
        HEAP8[((pdestEnd+i)|0)]=HEAP8[((psrc+i)|0)];
        i = (i+1)|0;
      } while (HEAP8[(((psrc)+(i-1))|0)]);
      return pdest|0;
    }
  function _strncmp(px, py, n) {
      var i = 0;
      while (i < n) {
        var x = HEAPU8[(((px)+(i))|0)];
        var y = HEAPU8[(((py)+(i))|0)];
        if (x == y && x == 0) return 0;
        if (x == 0) return -1;
        if (y == 0) return 1;
        if (x == y) {
          i ++;
          continue;
        } else {
          return x > y ? 1 : -1;
        }
      }
      return 0;
    }function _strcmp(px, py) {
      return _strncmp(px, py, TOTAL_MEMORY);
    }
  function _sprintf(s, format, varargs) {
      // int sprintf(char *restrict s, const char *restrict format, ...);
      // http://pubs.opengroup.org/onlinepubs/000095399/functions/printf.html
      return _snprintf(s, undefined, format, varargs);
    }
  function _fgetc(stream) {
      // int fgetc(FILE *stream);
      // http://pubs.opengroup.org/onlinepubs/000095399/functions/fgetc.html
      if (!FS.streams[stream]) return -1;
      var streamObj = FS.streams[stream];
      if (streamObj.eof || streamObj.error) return -1;
      var ret = _read(stream, _fgetc.ret, 1);
      if (ret == 0) {
        streamObj.eof = true;
        return -1;
      } else if (ret == -1) {
        streamObj.error = true;
        return -1;
      } else {
        return HEAPU8[((_fgetc.ret)|0)];
      }
    }function _fgets(s, n, stream) {
      // char *fgets(char *restrict s, int n, FILE *restrict stream);
      // http://pubs.opengroup.org/onlinepubs/000095399/functions/fgets.html
      if (!FS.streams[stream]) return 0;
      var streamObj = FS.streams[stream];
      if (streamObj.error || streamObj.eof) return 0;
      var byte_;
      for (var i = 0; i < n - 1 && byte_ != 10; i++) {
        byte_ = _fgetc(stream);
        if (byte_ == -1) {
          if (streamObj.error || (streamObj.eof && i == 0)) return 0;
          else if (streamObj.eof) break;
        }
        HEAP8[(((s)+(i))|0)]=byte_
      }
      HEAP8[(((s)+(i))|0)]=0
      return s;
    }
  function ___errno_location() {
      return ___errno_state;
    }var ___errno=___errno_location;
  function __isFloat(text) {
      return !!(/^[+-]?[0-9]*\.?[0-9]+([eE][+-]?[0-9]+)?$/.exec(text));
    }function __scanString(format, get, unget, varargs) {
      if (!__scanString.whiteSpace) {
        __scanString.whiteSpace = {};
        __scanString.whiteSpace[32] = 1;
        __scanString.whiteSpace[9] = 1;
        __scanString.whiteSpace[10] = 1;
        __scanString.whiteSpace[11] = 1;
        __scanString.whiteSpace[12] = 1;
        __scanString.whiteSpace[13] = 1;
        __scanString.whiteSpace[' '] = 1;
        __scanString.whiteSpace['\t'] = 1;
        __scanString.whiteSpace['\n'] = 1;
        __scanString.whiteSpace['\v'] = 1;
        __scanString.whiteSpace['\f'] = 1;
        __scanString.whiteSpace['\r'] = 1;
      }
      // Supports %x, %4x, %d.%d, %lld, %s, %f, %lf.
      // TODO: Support all format specifiers.
      format = Pointer_stringify(format);
      var soFar = 0;
      if (format.indexOf('%n') >= 0) {
        // need to track soFar
        var _get = get;
        get = function() {
          soFar++;
          return _get();
        }
        var _unget = unget;
        unget = function() {
          soFar--;
          return _unget();
        }
      }
      var formatIndex = 0;
      var argsi = 0;
      var fields = 0;
      var argIndex = 0;
      var next;
      mainLoop:
      for (var formatIndex = 0; formatIndex < format.length;) {
        if (format[formatIndex] === '%' && format[formatIndex+1] == 'n') {
          var argPtr = HEAP32[(((varargs)+(argIndex))>>2)];
          argIndex += Runtime.getAlignSize('void*', null, true);
          HEAP32[((argPtr)>>2)]=soFar;
          formatIndex += 2;
          continue;
        }
        // TODO: Support strings like "%5c" etc.
        if (format[formatIndex] === '%' && format[formatIndex+1] == 'c') {
          var argPtr = HEAP32[(((varargs)+(argIndex))>>2)];
          argIndex += Runtime.getAlignSize('void*', null, true);
          fields++;
          next = get();
          HEAP8[(argPtr)]=next
          formatIndex += 2;
          continue;
        }
        // remove whitespace
        while (1) {
          next = get();
          if (next == 0) return fields;
          if (!(next in __scanString.whiteSpace)) break;
        }
        unget();
        if (format[formatIndex] === '%') {
          formatIndex++;
          var suppressAssignment = false;
          if (format[formatIndex] == '*') {
            suppressAssignment = true;
            formatIndex++;
          }
          var maxSpecifierStart = formatIndex;
          while (format[formatIndex].charCodeAt(0) >= 48 &&
                 format[formatIndex].charCodeAt(0) <= 57) {
            formatIndex++;
          }
          var max_;
          if (formatIndex != maxSpecifierStart) {
            max_ = parseInt(format.slice(maxSpecifierStart, formatIndex), 10);
          }
          var long_ = false;
          var half = false;
          var longLong = false;
          if (format[formatIndex] == 'l') {
            long_ = true;
            formatIndex++;
            if (format[formatIndex] == 'l') {
              longLong = true;
              formatIndex++;
            }
          } else if (format[formatIndex] == 'h') {
            half = true;
            formatIndex++;
          }
          var type = format[formatIndex];
          formatIndex++;
          var curr = 0;
          var buffer = [];
          // Read characters according to the format. floats are trickier, they may be in an unfloat state in the middle, then be a valid float later
          if (type == 'f' || type == 'e' || type == 'g' ||
              type == 'F' || type == 'E' || type == 'G') {
            var last = 0;
            next = get();
            while (next > 0) {
              buffer.push(String.fromCharCode(next));
              if (__isFloat(buffer.join(''))) {
                last = buffer.length;
              }
              next = get();
            }
            for (var i = 0; i < buffer.length - last + 1; i++) {
              unget();
            }
            buffer.length = last;
          } else {
            next = get();
            var first = true;
            while ((curr < max_ || isNaN(max_)) && next > 0) {
              if (!(next in __scanString.whiteSpace) && // stop on whitespace
                  (type == 's' ||
                   ((type === 'd' || type == 'u' || type == 'i') && ((next >= 48 && next <= 57) ||
                                                                     (first && next == 45))) ||
                   ((type === 'x' || type === 'X') && (next >= 48 && next <= 57 ||
                                     next >= 97 && next <= 102 ||
                                     next >= 65 && next <= 70))) &&
                  (formatIndex >= format.length || next !== format[formatIndex].charCodeAt(0))) { // Stop when we read something that is coming up
                buffer.push(String.fromCharCode(next));
                next = get();
                curr++;
                first = false;
              } else {
                break;
              }
            }
            unget();
          }
          if (buffer.length === 0) return 0;  // Failure.
          if (suppressAssignment) continue;
          var text = buffer.join('');
          var argPtr = HEAP32[(((varargs)+(argIndex))>>2)];
          argIndex += Runtime.getAlignSize('void*', null, true);
          switch (type) {
            case 'd': case 'u': case 'i':
              if (half) {
                HEAP16[((argPtr)>>1)]=parseInt(text, 10);
              } else if (longLong) {
                (tempI64 = [parseInt(text, 10)>>>0,Math.min(Math.floor((parseInt(text, 10))/4294967296), 4294967295)>>>0],HEAP32[((argPtr)>>2)]=tempI64[0],HEAP32[(((argPtr)+(4))>>2)]=tempI64[1]);
              } else {
                HEAP32[((argPtr)>>2)]=parseInt(text, 10);
              }
              break;
            case 'X':
            case 'x':
              HEAP32[((argPtr)>>2)]=parseInt(text, 16)
              break;
            case 'F':
            case 'f':
            case 'E':
            case 'e':
            case 'G':
            case 'g':
            case 'E':
              // fallthrough intended
              if (long_) {
                HEAPF64[((argPtr)>>3)]=parseFloat(text)
              } else {
                HEAPF32[((argPtr)>>2)]=parseFloat(text)
              }
              break;
            case 's':
              var array = intArrayFromString(text);
              for (var j = 0; j < array.length; j++) {
                HEAP8[(((argPtr)+(j))|0)]=array[j]
              }
              break;
          }
          fields++;
        } else if (format[formatIndex] in __scanString.whiteSpace) {
          next = get();
          while (next in __scanString.whiteSpace) {
            if (next <= 0) break mainLoop;  // End of input.
            next = get();
          }
          unget(next);
          formatIndex++;
        } else {
          // Not a specifier.
          next = get();
          if (format[formatIndex].charCodeAt(0) !== next) {
            unget(next);
            break mainLoop;
          }
          formatIndex++;
        }
      }
      return fields;
    }function _sscanf(s, format, varargs) {
      // int sscanf(const char *restrict s, const char *restrict format, ... );
      // http://pubs.opengroup.org/onlinepubs/000095399/functions/scanf.html
      var index = 0;
      var get = function() { return HEAP8[(((s)+(index++))|0)]; };
      var unget = function() { index--; };
      return __scanString(format, get, unget, varargs);
    }
  function _strchr(ptr, chr) {
      ptr--;
      do {
        ptr++;
        var val = HEAP8[(ptr)];
        if (val == chr) return ptr;
      } while (val);
      return 0;
    }
  function _isalnum(chr) {
      return (chr >= 48 && chr <= 57) ||
             (chr >= 97 && chr <= 122) ||
             (chr >= 65 && chr <= 90);
    }
  function _isalpha(chr) {
      return (chr >= 97 && chr <= 122) ||
             (chr >= 65 && chr <= 90);
    }
  function _tolower(chr) {
      chr = chr|0;
      if ((chr|0) < 65) return chr|0;
      if ((chr|0) > 90) return chr|0;
      return (chr - 65 + 97)|0;
    }function _strncasecmp(px, py, n) {
      px = px|0; py = py|0; n = n|0;
      var i = 0, x = 0, y = 0;
      while ((i>>>0) < (n>>>0)) {
        x = _tolower(HEAP8[(((px)+(i))|0)])|0;
        y = _tolower(HEAP8[(((py)+(i))|0)])|0;
        if (((x|0) == (y|0)) & ((x|0) == 0)) return 0;
        if ((x|0) == 0) return -1;
        if ((y|0) == 0) return 1;
        if ((x|0) == (y|0)) {
          i = (i + 1)|0;
          continue;
        } else {
          return ((x>>>0) > (y>>>0) ? 1 : -1)|0;
        }
      }
      return 0;
    }function _strcasecmp(px, py) {
      px = px|0; py = py|0;
      return _strncasecmp(px, py, -1)|0;
    }
  var _environ=allocate(1, "i32*", ALLOC_STATIC);var ___environ=_environ;function ___buildEnvironment(env) {
      // WARNING: Arbitrary limit!
      var MAX_ENV_VALUES = 64;
      var TOTAL_ENV_SIZE = 1024;
      // Statically allocate memory for the environment.
      var poolPtr;
      var envPtr;
      if (!___buildEnvironment.called) {
        ___buildEnvironment.called = true;
        // Set default values. Use string keys for Closure Compiler compatibility.
        ENV['USER'] = 'root';
        ENV['PATH'] = '/';
        ENV['PWD'] = '/';
        ENV['HOME'] = '/home/emscripten';
        ENV['LANG'] = 'en_US.UTF-8';
        ENV['_'] = './this.program';
        // Allocate memory.
        poolPtr = allocate(TOTAL_ENV_SIZE, 'i8', ALLOC_STATIC);
        envPtr = allocate(MAX_ENV_VALUES * 4,
                          'i8*', ALLOC_STATIC);
        HEAP32[((envPtr)>>2)]=poolPtr
        HEAP32[((_environ)>>2)]=envPtr;
      } else {
        envPtr = HEAP32[((_environ)>>2)];
        poolPtr = HEAP32[((envPtr)>>2)];
      }
      // Collect key=value lines.
      var strings = [];
      var totalSize = 0;
      for (var key in env) {
        if (typeof env[key] === 'string') {
          var line = key + '=' + env[key];
          strings.push(line);
          totalSize += line.length;
        }
      }
      if (totalSize > TOTAL_ENV_SIZE) {
        throw new Error('Environment size exceeded TOTAL_ENV_SIZE!');
      }
      // Make new.
      var ptrSize = 4;
      for (var i = 0; i < strings.length; i++) {
        var line = strings[i];
        for (var j = 0; j < line.length; j++) {
          HEAP8[(((poolPtr)+(j))|0)]=line.charCodeAt(j);
        }
        HEAP8[(((poolPtr)+(j))|0)]=0;
        HEAP32[(((envPtr)+(i * ptrSize))>>2)]=poolPtr;
        poolPtr += line.length + 1;
      }
      HEAP32[(((envPtr)+(strings.length * ptrSize))>>2)]=0;
    }var ENV={};function _getenv(name) {
      // char *getenv(const char *name);
      // http://pubs.opengroup.org/onlinepubs/009695399/functions/getenv.html
      if (name === 0) return 0;
      name = Pointer_stringify(name);
      if (!ENV.hasOwnProperty(name)) return 0;
      if (_getenv.ret) _free(_getenv.ret);
      _getenv.ret = allocate(intArrayFromString(ENV[name]), 'i8', ALLOC_NORMAL);
      return _getenv.ret;
    }
  function _feof(stream) {
      // int feof(FILE *stream);
      // http://pubs.opengroup.org/onlinepubs/000095399/functions/feof.html
      return Number(FS.streams[stream] && FS.streams[stream].eof);
    }
  function _strstr(ptr1, ptr2) {
      var check = 0, start;
      do {
        if (!check) {
          start = ptr1;
          check = ptr2;
        }
        var curr1 = HEAP8[((ptr1++)|0)];
        var curr2 = HEAP8[((check++)|0)];
        if (curr2 == 0) return start;
        if (curr2 != curr1) {
          // rewind to one character after start, to find ez in eeez
          ptr1 = start + 1;
          check = 0;
        }
      } while (curr1);
      return 0;
    }
  function _strrchr(ptr, chr) {
      var ptr2 = ptr + _strlen(ptr);
      do {
        if (HEAP8[(ptr2)] == chr) return ptr2;
        ptr2--;
      } while (ptr2 >= ptr);
      return 0;
    }
  var ___stat_struct_layout={__size__:68,st_dev:0,st_ino:4,st_mode:8,st_nlink:12,st_uid:16,st_gid:20,st_rdev:24,st_size:28,st_atime:32,st_spare1:36,st_mtime:40,st_spare2:44,st_ctime:48,st_spare3:52,st_blksize:56,st_blocks:60,st_spare4:64};function _stat(path, buf, dontResolveLastLink) {
      // http://pubs.opengroup.org/onlinepubs/7908799/xsh/stat.html
      // int stat(const char *path, struct stat *buf);
      // NOTE: dontResolveLastLink is a shortcut for lstat(). It should never be
      //       used in client code.
      var obj = FS.findObject(Pointer_stringify(path), dontResolveLastLink);
      if (obj === null || !FS.forceLoadFile(obj)) return -1;
      var offsets = ___stat_struct_layout;
      // Constants.
      HEAP32[(((buf)+(offsets.st_nlink))>>2)]=1
      HEAP32[(((buf)+(offsets.st_uid))>>2)]=0
      HEAP32[(((buf)+(offsets.st_gid))>>2)]=0
      HEAP32[(((buf)+(offsets.st_blksize))>>2)]=4096
      // Variables.
      HEAP32[(((buf)+(offsets.st_ino))>>2)]=obj.inodeNumber
      var time = Math.floor(obj.timestamp / 1000);
      if (offsets.st_atime === undefined) {
        offsets.st_atime = offsets.st_atim.tv_sec;
        offsets.st_mtime = offsets.st_mtim.tv_sec;
        offsets.st_ctime = offsets.st_ctim.tv_sec;
        var nanosec = (obj.timestamp % 1000) * 1000;
        HEAP32[(((buf)+(offsets.st_atim.tv_nsec))>>2)]=nanosec
        HEAP32[(((buf)+(offsets.st_mtim.tv_nsec))>>2)]=nanosec
        HEAP32[(((buf)+(offsets.st_ctim.tv_nsec))>>2)]=nanosec
      }
      HEAP32[(((buf)+(offsets.st_atime))>>2)]=time
      HEAP32[(((buf)+(offsets.st_mtime))>>2)]=time
      HEAP32[(((buf)+(offsets.st_ctime))>>2)]=time
      var mode = 0;
      var size = 0;
      var blocks = 0;
      var dev = 0;
      var rdev = 0;
      if (obj.isDevice) {
        //  Device numbers reuse inode numbers.
        dev = rdev = obj.inodeNumber;
        size = blocks = 0;
        mode = 0x2000;  // S_IFCHR.
      } else {
        dev = 1;
        rdev = 0;
        // NOTE: In our implementation, st_blocks = Math.ceil(st_size/st_blksize),
        //       but this is not required by the standard.
        if (obj.isFolder) {
          size = 4096;
          blocks = 1;
          mode = 0x4000;  // S_IFDIR.
        } else {
          var data = obj.contents || obj.link;
          size = data.length;
          blocks = Math.ceil(data.length / 4096);
          mode = obj.link === undefined ? 0x8000 : 0xA000;  // S_IFREG, S_IFLNK.
        }
      }
      HEAP32[(((buf)+(offsets.st_dev))>>2)]=dev;
      HEAP32[(((buf)+(offsets.st_rdev))>>2)]=rdev;
      HEAP32[(((buf)+(offsets.st_size))>>2)]=size
      HEAP32[(((buf)+(offsets.st_blocks))>>2)]=blocks
      if (obj.read) mode |= 0x16D;  // S_IRUSR | S_IXUSR | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH.
      if (obj.write) mode |= 0x92;  // S_IWUSR | S_IWGRP | S_IWOTH.
      HEAP32[(((buf)+(offsets.st_mode))>>2)]=mode
      return 0;
    }
  function _isspace(chr) {
      return chr in { 32: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0 };
    }function __parseInt(str, endptr, base, min, max, bits, unsign) {
      // Skip space.
      while (_isspace(HEAP8[(str)])) str++;
      // Check for a plus/minus sign.
      var multiplier = 1;
      if (HEAP8[(str)] == 45) {
        multiplier = -1;
        str++;
      } else if (HEAP8[(str)] == 43) {
        str++;
      }
      // Find base.
      var finalBase = base;
      if (!finalBase) {
        if (HEAP8[(str)] == 48) {
          if (HEAP8[((str+1)|0)] == 120 ||
              HEAP8[((str+1)|0)] == 88) {
            finalBase = 16;
            str += 2;
          } else {
            finalBase = 8;
            str++;
          }
        }
      } else if (finalBase==16) {
        if (HEAP8[(str)] == 48) {
          if (HEAP8[((str+1)|0)] == 120 ||
              HEAP8[((str+1)|0)] == 88) {
            str += 2;
          }
        }
      }
      if (!finalBase) finalBase = 10;
      // Get digits.
      var chr;
      var ret = 0;
      while ((chr = HEAP8[(str)]) != 0) {
        var digit = parseInt(String.fromCharCode(chr), finalBase);
        if (isNaN(digit)) {
          break;
        } else {
          ret = ret * finalBase + digit;
          str++;
        }
      }
      // Apply sign.
      ret *= multiplier;
      // Set end pointer.
      if (endptr) {
        HEAP32[((endptr)>>2)]=str
      }
      // Unsign if needed.
      if (unsign) {
        if (Math.abs(ret) > max) {
          ret = max;
          ___setErrNo(ERRNO_CODES.ERANGE);
        } else {
          ret = unSign(ret, bits);
        }
      }
      // Validate range.
      if (ret > max || ret < min) {
        ret = ret > max ? max : min;
        ___setErrNo(ERRNO_CODES.ERANGE);
      }
      if (bits == 64) {
        return tempRet0 = Math.min(Math.floor((ret)/4294967296), 4294967295)>>>0,ret>>>0;
      }
      return ret;
    }function _strtol(str, endptr, base) {
      return __parseInt(str, endptr, base, -2147483648, 2147483647, 32);  // LONG_MIN, LONG_MAX.
    }function _atoi(ptr) {
      return _strtol(ptr, null, 10);
    }
;
;
  function _strdup(ptr) {
      var len = _strlen(ptr);
      var newStr = _malloc(len + 1);
      (_memcpy(newStr, ptr, len)|0);
      HEAP8[(((newStr)+(len))|0)]=0;
      return newStr;
    }
  var ERRNO_MESSAGES={0:"Success",1:"Not super-user",2:"No such file or directory",3:"No such process",4:"Interrupted system call",5:"I/O error",6:"No such device or address",7:"Arg list too long",8:"Exec format error",9:"Bad file number",10:"No children",11:"No more processes",12:"Not enough core",13:"Permission denied",14:"Bad address",15:"Block device required",16:"Mount device busy",17:"File exists",18:"Cross-device link",19:"No such device",20:"Not a directory",21:"Is a directory",22:"Invalid argument",23:"Too many open files in system",24:"Too many open files",25:"Not a typewriter",26:"Text file busy",27:"File too large",28:"No space left on device",29:"Illegal seek",30:"Read only file system",31:"Too many links",32:"Broken pipe",33:"Math arg out of domain of func",34:"Math result not representable",35:"No message of desired type",36:"Identifier removed",37:"Channel number out of range",38:"Level 2 not synchronized",39:"Level 3 halted",40:"Level 3 reset",41:"Link number out of range",42:"Protocol driver not attached",43:"No CSI structure available",44:"Level 2 halted",45:"Deadlock condition",46:"No record locks available",50:"Invalid exchange",51:"Invalid request descriptor",52:"Exchange full",53:"No anode",54:"Invalid request code",55:"Invalid slot",56:"File locking deadlock error",57:"Bad font file fmt",60:"Device not a stream",61:"No data (for no delay io)",62:"Timer expired",63:"Out of streams resources",64:"Machine is not on the network",65:"Package not installed",66:"The object is remote",67:"The link has been severed",68:"Advertise error",69:"Srmount error",70:"Communication error on send",71:"Protocol error",74:"Multihop attempted",75:"Inode is remote (not really error)",76:"Cross mount point (not really error)",77:"Trying to read unreadable message",79:"Inappropriate file type or format",80:"Given log. name not unique",81:"f.d. invalid for this operation",82:"Remote address changed",83:"Can\t access a needed shared lib",84:"Accessing a corrupted shared lib",85:".lib section in a.out corrupted",86:"Attempting to link in too many libs",87:"Attempting to exec a shared library",88:"Function not implemented",89:"No more files",90:"Directory not empty",91:"File or path name too long",92:"Too many symbolic links",95:"Operation not supported on transport endpoint",96:"Protocol family not supported",104:"Connection reset by peer",105:"No buffer space available",106:"Address family not supported by protocol family",107:"Protocol wrong type for socket",108:"Socket operation on non-socket",109:"Protocol not available",110:"Can't send after socket shutdown",111:"Connection refused",112:"Address already in use",113:"Connection aborted",114:"Network is unreachable",115:"Network interface is not configured",116:"Connection timed out",117:"Host is down",118:"Host is unreachable",119:"Connection already in progress",120:"Socket already connected",121:"Destination address required",122:"Message too long",123:"Unknown protocol",124:"Socket type not supported",125:"Address not available",126:"ENETRESET",127:"Socket is already connected",128:"Socket is not connected",129:"TOOMANYREFS",130:"EPROCLIM",131:"EUSERS",132:"EDQUOT",133:"ESTALE",134:"Not supported",135:"No medium (in tape drive)",136:"No such host or network path",137:"Filename exists with different case",138:"EILSEQ",139:"Value too large for defined data type",140:"Operation canceled",141:"State not recoverable",142:"Previous owner died",143:"Streams pipe error"};function _strerror_r(errnum, strerrbuf, buflen) {
      if (errnum in ERRNO_MESSAGES) {
        if (ERRNO_MESSAGES[errnum].length > buflen - 1) {
          return ___setErrNo(ERRNO_CODES.ERANGE);
        } else {
          var msg = ERRNO_MESSAGES[errnum];
          for (var i = 0; i < msg.length; i++) {
            HEAP8[(((strerrbuf)+(i))|0)]=msg.charCodeAt(i)
          }
          HEAP8[(((strerrbuf)+(i))|0)]=0
          return 0;
        }
      } else {
        return ___setErrNo(ERRNO_CODES.EINVAL);
      }
    }function _strerror(errnum) {
      if (!_strerror.buffer) _strerror.buffer = _malloc(256);
      _strerror_r(errnum, _strerror.buffer, 256);
      return _strerror.buffer;
    }
  function _vsprintf(s, format, va_arg) {
      return _sprintf(s, format, HEAP32[((va_arg)>>2)]);
    }
  function _strncpy(pdest, psrc, num) {
      pdest = pdest|0; psrc = psrc|0; num = num|0;
      var padding = 0, curr = 0, i = 0;
      while ((i|0) < (num|0)) {
        curr = padding ? 0 : HEAP8[(((psrc)+(i))|0)];
        HEAP8[(((pdest)+(i))|0)]=curr
        padding = padding ? 1 : (HEAP8[(((psrc)+(i))|0)] == 0);
        i = (i+1)|0;
      }
      return pdest|0;
    }
;
;
;
;
  function _bsearch(key, base, num, size, compar) {
      var cmp = function(x, y) {
        return Runtime.dynCall('iii', compar, [x, y])
      };
      var left = 0;
      var right = num;
      var mid, test, addr;
      while (left < right) {
        mid = (left + right) >>> 1;
        addr = base + (mid * size);
        test = cmp(key, addr);
        if (test < 0) {
          right = mid;
        } else if (test > 0) {
          left = mid + 1;
        } else {
          return addr;
        }
      }
      return 0;
    }
  function _memcmp(p1, p2, num) {
      p1 = p1|0; p2 = p2|0; num = num|0;
      var i = 0, v1 = 0, v2 = 0;
      while ((i|0) < (num|0)) {
        var v1 = HEAPU8[(((p1)+(i))|0)];
        var v2 = HEAPU8[(((p2)+(i))|0)];
        if ((v1|0) != (v2|0)) return ((v1|0) > (v2|0) ? 1 : -1)|0;
        i = (i+1)|0;
      }
      return 0;
    }
  var _sqrt=Math.sqrt;
  function _isupper(chr) {
      return chr >= 65 && chr <= 90;
    }
  var _cos=Math.cos;
  var _sin=Math.sin;
  var _tan=Math.tan;
  var _exp=Math.exp;
  var _atan2=Math.atan2;
  var _floor=Math.floor;
  var _atan=Math.atan;
;
  function _setlocale(category, locale) {
      if (!_setlocale.ret) _setlocale.ret = allocate([0], 'i8', ALLOC_NORMAL);
      return _setlocale.ret;
    }
  function _strtok_r(s, delim, lasts) {
      var skip_leading_delim = 1;
      var spanp;
      var c, sc;
      var tok;
      if (s == 0 && (s = getValue(lasts, 'i8*')) == 0) {
        return 0;
      }
      cont: while (1) {
        c = getValue(s++, 'i8');
        for (spanp = delim; (sc = getValue(spanp++, 'i8')) != 0;) {
          if (c == sc) {
            if (skip_leading_delim) {
              continue cont;
            } else {
              setValue(lasts, s, 'i8*');
              setValue(s - 1, 0, 'i8');
              return s - 1;
            }
          }
        }
        break;
      }
      if (c == 0) {
        setValue(lasts, 0, 'i8*');
        return 0;
      }
      tok = s - 1;
      for (;;) {
        c = getValue(s++, 'i8');
        spanp = delim;
        do {
          if ((sc = getValue(spanp++, 'i8')) == c) {
            if (c == 0) {
              s = 0;
            } else {
              setValue(s - 1, 0, 'i8');
            }
            setValue(lasts, s, 'i8*');
            return tok;
          }
        } while (sc != 0);
      }
      abort('strtok_r error!');
    }
  function _strpbrk(ptr1, ptr2) {
      var curr;
      var searchSet = {};
      while (1) {
        var curr = HEAP8[((ptr2++)|0)];
        if (!curr) break;
        searchSet[curr] = 1;
      }
      while (1) {
        curr = HEAP8[(ptr1)];
        if (!curr) break;
        if (curr in searchSet) return ptr1;
        ptr1++;
      }
      return 0;
    }
  var ___strtok_state=0;function _strtok(s, delim) {
      return _strtok_r(s, delim, ___strtok_state);
    }
  var _fabs=Math.abs;
  function _putenv(string) {
      // int putenv(char *string);
      // http://pubs.opengroup.org/onlinepubs/009695399/functions/putenv.html
      // WARNING: According to the standard (and the glibc implementation), the
      //          string is taken by reference so future changes are reflected.
      //          We copy it instead, possibly breaking some uses.
      if (string === 0) {
        ___setErrNo(ERRNO_CODES.EINVAL);
        return -1;
      }
      string = Pointer_stringify(string);
      var splitPoint = string.indexOf('=')
      if (string === '' || string.indexOf('=') === -1) {
        ___setErrNo(ERRNO_CODES.EINVAL);
        return -1;
      }
      var name = string.slice(0, splitPoint);
      var value = string.slice(splitPoint + 1);
      if (!(name in ENV) || ENV[name] !== value) {
        ENV[name] = value;
        ___buildEnvironment(ENV);
      }
      return 0;
    }
;
  var _setjmp=undefined;
  function _longjmp(env, value) {
      throw { longjmp: true, id: HEAP32[((env)>>2)], value: value || 1 };
    }
;
  function _fstat(fildes, buf) {
      // int fstat(int fildes, struct stat *buf);
      // http://pubs.opengroup.org/onlinepubs/7908799/xsh/fstat.html
      if (!FS.streams[fildes]) {
        ___setErrNo(ERRNO_CODES.EBADF);
        return -1;
      } else {
        var pathArray = intArrayFromString(FS.streams[fildes].path);
        return _stat(allocate(pathArray, 'i8', ALLOC_STACK), buf);
      }
    }
  function _fileno(stream) {
      // int fileno(FILE *stream);
      // http://pubs.opengroup.org/onlinepubs/000095399/functions/fileno.html
      // We use file descriptor numbers and FILE* streams interchangeably.
      return stream;
    }
  function _hypot(a, b) {
       return Math.sqrt(a*a + b*b);
    }
  function _rand() {
      return Math.floor(Math.random()*0x80000000);
    }
  function _access(path, amode) {
      // int access(const char *path, int amode);
      // http://pubs.opengroup.org/onlinepubs/000095399/functions/access.html
      path = Pointer_stringify(path);
      var target = FS.findObject(path);
      if (target === null) return -1;
      if ((amode & 2 && !target.write) ||  // W_OK.
          ((amode & 1 || amode & 4) && !target.read)) {  // X_OK, R_OK.
        ___setErrNo(ERRNO_CODES.EACCES);
        return -1;
      } else {
        return 0;
      }
    }
  function _cbrt(x) {
      return Math.pow(x, 1/3);
    }
  var _asin=Math.asin;
;
;
;
  function _qsort(base, num, size, cmp) {
      if (num == 0 || size == 0) return;
      // forward calls to the JavaScript sort method
      // first, sort the items logically
      var comparator = function(x, y) {
        return Runtime.dynCall('iii', cmp, [x, y]);
      }
      var keys = [];
      for (var i = 0; i < num; i++) keys.push(i);
      keys.sort(function(a, b) {
        return comparator(base+a*size, base+b*size);
      });
      // apply the sort
      var temp = _malloc(num*size);
      _memcpy(temp, base, num*size);
      for (var i = 0; i < num; i++) {
        if (keys[i] == i) continue; // already in place
        _memcpy(base+i*size, temp+keys[i]*size, size);
      }
      _free(temp);
    }
;
  var _ceil=Math.ceil;
;
;
  function _getgid() {
      // gid_t getgid(void);
      // http://pubs.opengroup.org/onlinepubs/000095399/functions/getgid.html
      // We have just one process/group/user, all with ID 0.
      return 0;
    }var _getpid=_getgid;
  function _time(ptr) {
      var ret = Math.floor(Date.now()/1000);
      if (ptr) {
        HEAP32[((ptr)>>2)]=ret
      }
      return ret;
    }
  function _srand(seed) {}
;
;
;
;
;
;
;
;
;
;
;
;
;
;
;
;
;
;
;
;
;
;
  var _llvm_pow_f64=Math.pow;
;
;
;
;
;
;
;
;
;
;
;
;
  function _mmap(start, num, prot, flags, stream, offset) {
      /* FIXME: Since mmap is normally implemented at the kernel level,
       * this implementation simply uses malloc underneath the call to
       * mmap.
       */
      var MAP_PRIVATE = 2;
      var allocated = false;
      if (!_mmap.mappings) _mmap.mappings = {};
      if (stream == -1) {
        var ptr = _malloc(num);
        if (!ptr) return -1;
        _memset(ptr, 0, num);
        allocated = true;
      } else {
        var info = FS.streams[stream];
        if (!info) return -1;
        var contents = info.object.contents;
        // Only make a new copy when MAP_PRIVATE is specified.
        if (flags & MAP_PRIVATE == 0) {
          // We can't emulate MAP_SHARED when the file is not backed by HEAP.
          assert(contents.buffer === HEAPU8.buffer);
          ptr = contents.byteOffset;
          allocated = false;
        } else {
          // Try to avoid unnecessary slices.
          if (offset > 0 || offset + num < contents.length) {
            if (contents.subarray) {
              contents = contents.subarray(offset, offset+num);
            } else {
              contents = Array.prototype.slice.call(contents, offset, offset+num);
            }
          }
          ptr = _malloc(num);
          if (!ptr) return -1;
          HEAPU8.set(contents, ptr);
          allocated = true;
        }
      }
      _mmap.mappings[ptr] = { malloc: ptr, num: num, allocated: allocated };
      return ptr;
    }
  function _munmap(start, num) {
      if (!_mmap.mappings) _mmap.mappings = {};
      // TODO: support unmmap'ing parts of allocations
      var info = _mmap.mappings[start];
      if (!info) return 0;
      if (num == info.num) {
        _mmap.mappings[start] = null;
        if (info.allocated) {
          _free(info.malloc);
        }
      }
      return 0;
    }
  function _log10(x) {
      return Math.log(x) / Math.LN10;
    }
  var _atanf=Math.atan;
;
;
;
;
;
;
;
;
  function _pthread_mutex_lock() {}
  function _pthread_mutex_unlock() {}
  function _pthread_cond_broadcast() {
      return 0;
    }
  function _pthread_cond_wait() {
      return 0;
    }
  function _atexit(func, arg) {
      __ATEXIT__.unshift({ func: func, arg: arg });
    }var ___cxa_atexit=_atexit;
  function ___cxa_allocate_exception(size) {
      return _malloc(size);
    }
  function ___cxa_free_exception(ptr) {
      try {
        return _free(ptr);
      } catch(e) { // XXX FIXME
      }
    }
  function _llvm_eh_exception() {
      return HEAP32[((_llvm_eh_exception.buf)>>2)];
    }
  function __ZSt18uncaught_exceptionv() { // std::uncaught_exception()
      return !!__ZSt18uncaught_exceptionv.uncaught_exception;
    }
  function ___cxa_is_number_type(type) {
      var isNumber = false;
      try { if (type == __ZTIi) isNumber = true } catch(e){}
      try { if (type == __ZTIj) isNumber = true } catch(e){}
      try { if (type == __ZTIl) isNumber = true } catch(e){}
      try { if (type == __ZTIm) isNumber = true } catch(e){}
      try { if (type == __ZTIx) isNumber = true } catch(e){}
      try { if (type == __ZTIy) isNumber = true } catch(e){}
      try { if (type == __ZTIf) isNumber = true } catch(e){}
      try { if (type == __ZTId) isNumber = true } catch(e){}
      try { if (type == __ZTIe) isNumber = true } catch(e){}
      try { if (type == __ZTIc) isNumber = true } catch(e){}
      try { if (type == __ZTIa) isNumber = true } catch(e){}
      try { if (type == __ZTIh) isNumber = true } catch(e){}
      try { if (type == __ZTIs) isNumber = true } catch(e){}
      try { if (type == __ZTIt) isNumber = true } catch(e){}
      return isNumber;
    }function ___cxa_does_inherit(definiteType, possibilityType, possibility) {
      if (possibility == 0) return false;
      if (possibilityType == 0 || possibilityType == definiteType)
        return true;
      var possibility_type_info;
      if (___cxa_is_number_type(possibilityType)) {
        possibility_type_info = possibilityType;
      } else {
        var possibility_type_infoAddr = HEAP32[((possibilityType)>>2)] - 8;
        possibility_type_info = HEAP32[((possibility_type_infoAddr)>>2)];
      }
      switch (possibility_type_info) {
      case 0: // possibility is a pointer
        // See if definite type is a pointer
        var definite_type_infoAddr = HEAP32[((definiteType)>>2)] - 8;
        var definite_type_info = HEAP32[((definite_type_infoAddr)>>2)];
        if (definite_type_info == 0) {
          // Also a pointer; compare base types of pointers
          var defPointerBaseAddr = definiteType+8;
          var defPointerBaseType = HEAP32[((defPointerBaseAddr)>>2)];
          var possPointerBaseAddr = possibilityType+8;
          var possPointerBaseType = HEAP32[((possPointerBaseAddr)>>2)];
          return ___cxa_does_inherit(defPointerBaseType, possPointerBaseType, possibility);
        } else
          return false; // one pointer and one non-pointer
      case 1: // class with no base class
        return false;
      case 2: // class with base class
        var parentTypeAddr = possibilityType + 8;
        var parentType = HEAP32[((parentTypeAddr)>>2)];
        return ___cxa_does_inherit(definiteType, parentType, possibility);
      default:
        return false; // some unencountered type
      }
    }
  function ___resumeException(ptr) {
      if (HEAP32[((_llvm_eh_exception.buf)>>2)] == 0) HEAP32[((_llvm_eh_exception.buf)>>2)]=ptr;
      throw ptr + " - Exception catching is disabled, this exception cannot be caught. Compile with -s DISABLE_EXCEPTION_CATCHING=0 or DISABLE_EXCEPTION_CATCHING=2 to catch.";;
    }function ___cxa_find_matching_catch(thrown, throwntype) {
      if (thrown == -1) thrown = HEAP32[((_llvm_eh_exception.buf)>>2)];
      if (throwntype == -1) throwntype = HEAP32[(((_llvm_eh_exception.buf)+(4))>>2)];
      var typeArray = Array.prototype.slice.call(arguments, 2);
      // If throwntype is a pointer, this means a pointer has been
      // thrown. When a pointer is thrown, actually what's thrown
      // is a pointer to the pointer. We'll dereference it.
      if (throwntype != 0 && !___cxa_is_number_type(throwntype)) {
        var throwntypeInfoAddr= HEAP32[((throwntype)>>2)] - 8;
        var throwntypeInfo= HEAP32[((throwntypeInfoAddr)>>2)];
        if (throwntypeInfo == 0)
          thrown = HEAP32[((thrown)>>2)];
      }
      // The different catch blocks are denoted by different types.
      // Due to inheritance, those types may not precisely match the
      // type of the thrown object. Find one which matches, and
      // return the type of the catch block which should be called.
      for (var i = 0; i < typeArray.length; i++) {
        if (___cxa_does_inherit(typeArray[i], throwntype, thrown))
          return tempRet0 = typeArray[i],thrown;
      }
      // Shouldn't happen unless we have bogus data in typeArray
      // or encounter a type for which emscripten doesn't have suitable
      // typeinfo defined. Best-efforts match just in case.
      return tempRet0 = throwntype,thrown;
    }function ___cxa_throw(ptr, type, destructor) {
      if (!___cxa_throw.initialized) {
        try {
          HEAP32[((__ZTVN10__cxxabiv119__pointer_type_infoE)>>2)]=0; // Workaround for libcxxabi integration bug
        } catch(e){}
        try {
          HEAP32[((__ZTVN10__cxxabiv117__class_type_infoE)>>2)]=1; // Workaround for libcxxabi integration bug
        } catch(e){}
        try {
          HEAP32[((__ZTVN10__cxxabiv120__si_class_type_infoE)>>2)]=2; // Workaround for libcxxabi integration bug
        } catch(e){}
        ___cxa_throw.initialized = true;
      }
      HEAP32[((_llvm_eh_exception.buf)>>2)]=ptr
      HEAP32[(((_llvm_eh_exception.buf)+(4))>>2)]=type
      HEAP32[(((_llvm_eh_exception.buf)+(8))>>2)]=destructor
      if (!("uncaught_exception" in __ZSt18uncaught_exceptionv)) {
        __ZSt18uncaught_exceptionv.uncaught_exception = 1;
      } else {
        __ZSt18uncaught_exceptionv.uncaught_exception++;
      }
      throw ptr + " - Exception catching is disabled, this exception cannot be caught. Compile with -s DISABLE_EXCEPTION_CATCHING=0 or DISABLE_EXCEPTION_CATCHING=2 to catch.";;
    }
  function ___cxa_begin_catch(ptr) {
      __ZSt18uncaught_exceptionv.uncaught_exception--;
      return ptr;
    }
  function ___cxa_end_catch() {
      if (___cxa_end_catch.rethrown) {
        ___cxa_end_catch.rethrown = false;
        return;
      }
      // Clear state flag.
      __THREW__ = 0;
      // Clear type.
      HEAP32[(((_llvm_eh_exception.buf)+(4))>>2)]=0
      // Call destructor if one is registered then clear it.
      var ptr = HEAP32[((_llvm_eh_exception.buf)>>2)];
      var destructor = HEAP32[(((_llvm_eh_exception.buf)+(8))>>2)];
      if (destructor) {
        Runtime.dynCall('vi', destructor, [ptr]);
        HEAP32[(((_llvm_eh_exception.buf)+(8))>>2)]=0
      }
      // Free ptr if it isn't null.
      if (ptr) {
        ___cxa_free_exception(ptr);
        HEAP32[((_llvm_eh_exception.buf)>>2)]=0
      }
    }
  function _llvm_trap() {
      throw 'trap! ' + new Error().stack;
    }
  function _ungetc(c, stream) {
      // int ungetc(int c, FILE *stream);
      // http://pubs.opengroup.org/onlinepubs/000095399/functions/ungetc.html
      if (FS.streams[stream]) {
        c = unSign(c & 0xFF);
        FS.streams[stream].ungotten.push(c);
        return c;
      } else {
        return -1;
      }
    }
  var _getc=_fgetc;
  function _abort() {
      ABORT = true;
      throw 'abort() at ' + (new Error().stack);
    }
  function _memmove(dest, src, num) {
      dest = dest|0; src = src|0; num = num|0;
      if (((src|0) < (dest|0)) & ((dest|0) < ((src + num)|0))) {
        // Unlikely case: Copy backwards in a safe manner
        src = (src + num)|0;
        dest = (dest + num)|0;
        while ((num|0) > 0) {
          dest = (dest - 1)|0;
          src = (src - 1)|0;
          num = (num - 1)|0;
          HEAP8[(dest)]=HEAP8[(src)];
        }
      } else {
        _memcpy(dest, src, num) | 0;
      }
    }var _llvm_memmove_p0i8_p0i8_i32=_memmove;
  function ___cxa_rethrow() {
      ___cxa_end_catch.rethrown = true;
      throw HEAP32[((_llvm_eh_exception.buf)>>2)] + " - Exception catching is disabled, this exception cannot be caught. Compile with -s DISABLE_EXCEPTION_CATCHING=0 or DISABLE_EXCEPTION_CATCHING=2 to catch.";;
    }
  function _sysconf(name) {
      // long sysconf(int name);
      // http://pubs.opengroup.org/onlinepubs/009695399/functions/sysconf.html
      switch(name) {
        case 8: return PAGE_SIZE;
        case 54:
        case 56:
        case 21:
        case 61:
        case 63:
        case 22:
        case 67:
        case 23:
        case 24:
        case 25:
        case 26:
        case 27:
        case 69:
        case 28:
        case 101:
        case 70:
        case 71:
        case 29:
        case 30:
        case 199:
        case 75:
        case 76:
        case 32:
        case 43:
        case 44:
        case 80:
        case 46:
        case 47:
        case 45:
        case 48:
        case 49:
        case 42:
        case 82:
        case 33:
        case 7:
        case 108:
        case 109:
        case 107:
        case 112:
        case 119:
        case 121:
          return 200809;
        case 13:
        case 104:
        case 94:
        case 95:
        case 34:
        case 35:
        case 77:
        case 81:
        case 83:
        case 84:
        case 85:
        case 86:
        case 87:
        case 88:
        case 89:
        case 90:
        case 91:
        case 94:
        case 95:
        case 110:
        case 111:
        case 113:
        case 114:
        case 115:
        case 116:
        case 117:
        case 118:
        case 120:
        case 40:
        case 16:
        case 79:
        case 19:
          return -1;
        case 92:
        case 93:
        case 5:
        case 72:
        case 6:
        case 74:
        case 92:
        case 93:
        case 96:
        case 97:
        case 98:
        case 99:
        case 102:
        case 103:
        case 105:
          return 1;
        case 38:
        case 66:
        case 50:
        case 51:
        case 4:
          return 1024;
        case 15:
        case 64:
        case 41:
          return 32;
        case 55:
        case 37:
        case 17:
          return 2147483647;
        case 18:
        case 1:
          return 47839;
        case 59:
        case 57:
          return 99;
        case 68:
        case 58:
          return 2048;
        case 0: return 2097152;
        case 3: return 65536;
        case 14: return 32768;
        case 73: return 32767;
        case 39: return 16384;
        case 60: return 1000;
        case 106: return 700;
        case 52: return 256;
        case 62: return 255;
        case 2: return 100;
        case 65: return 64;
        case 36: return 20;
        case 100: return 16;
        case 20: return 6;
        case 53: return 4;
        case 10: return 1;
      }
      ___setErrNo(ERRNO_CODES.EINVAL);
      return -1;
    }
  function _isxdigit(chr) {
      return (chr >= 48 && chr <= 57) ||
             (chr >= 97 && chr <= 102) ||
             (chr >= 65 && chr <= 70);
    }var _isxdigit_l=_isxdigit;
  function _isdigit(chr) {
      return chr >= 48 && chr <= 57;
    }var _isdigit_l=_isdigit;
  function __Z7catopenPKci() { throw 'catopen not implemented' }
  function __Z7catgetsP8_nl_catdiiPKc() { throw 'catgets not implemented' }
  function __Z8catcloseP8_nl_catd() { throw 'catclose not implemented' }
  function _newlocale(mask, locale, base) {
      return 0;
    }
  function _freelocale(locale) {}
  function ___ctype_b_loc() {
      // http://refspecs.freestandards.org/LSB_3.0.0/LSB-Core-generic/LSB-Core-generic/baselib---ctype-b-loc.html
      var me = ___ctype_b_loc;
      if (!me.ret) {
        var values = [
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,2,2,2,2,2,2,2,2,2,8195,8194,8194,8194,8194,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,24577,49156,49156,49156,
          49156,49156,49156,49156,49156,49156,49156,49156,49156,49156,49156,49156,55304,55304,55304,55304,55304,55304,55304,55304,
          55304,55304,49156,49156,49156,49156,49156,49156,49156,54536,54536,54536,54536,54536,54536,50440,50440,50440,50440,50440,
          50440,50440,50440,50440,50440,50440,50440,50440,50440,50440,50440,50440,50440,50440,50440,49156,49156,49156,49156,49156,
          49156,54792,54792,54792,54792,54792,54792,50696,50696,50696,50696,50696,50696,50696,50696,50696,50696,50696,50696,50696,
          50696,50696,50696,50696,50696,50696,50696,49156,49156,49156,49156,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
        ];
        var i16size = 2;
        var arr = _malloc(values.length * i16size);
        for (var i = 0; i < values.length; i++) {
          HEAP16[(((arr)+(i * i16size))>>1)]=values[i]
        }
        me.ret = allocate([arr + 128 * i16size], 'i16*', ALLOC_NORMAL);
      }
      return me.ret;
    }
  function ___ctype_tolower_loc() {
      // http://refspecs.freestandards.org/LSB_3.1.1/LSB-Core-generic/LSB-Core-generic/libutil---ctype-tolower-loc.html
      var me = ___ctype_tolower_loc;
      if (!me.ret) {
        var values = [
          128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,
          158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,
          188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,
          218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,
          248,249,250,251,252,253,254,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,
          33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,97,98,99,100,101,102,103,
          104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,91,92,93,94,95,96,97,98,99,100,101,102,103,
          104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,
          134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,
          164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,
          194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,
          224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,
          254,255
        ];
        var i32size = 4;
        var arr = _malloc(values.length * i32size);
        for (var i = 0; i < values.length; i++) {
          HEAP32[(((arr)+(i * i32size))>>2)]=values[i]
        }
        me.ret = allocate([arr + 128 * i32size], 'i32*', ALLOC_NORMAL);
      }
      return me.ret;
    }
  function ___ctype_toupper_loc() {
      // http://refspecs.freestandards.org/LSB_3.1.1/LSB-Core-generic/LSB-Core-generic/libutil---ctype-toupper-loc.html
      var me = ___ctype_toupper_loc;
      if (!me.ret) {
        var values = [
          128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,
          158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,
          188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,
          218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,
          248,249,250,251,252,253,254,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,
          33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,
          73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,
          81,82,83,84,85,86,87,88,89,90,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,
          145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,
          175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,
          205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,
          235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255
        ];
        var i32size = 4;
        var arr = _malloc(values.length * i32size);
        for (var i = 0; i < values.length; i++) {
          HEAP32[(((arr)+(i * i32size))>>2)]=values[i]
        }
        me.ret = allocate([arr + 128 * i32size], 'i32*', ALLOC_NORMAL);
      }
      return me.ret;
    }
  function _strftime(s, maxsize, format, timeptr) {
      // size_t strftime(char *restrict s, size_t maxsize, const char *restrict format, const struct tm *restrict timeptr);
      // http://pubs.opengroup.org/onlinepubs/009695399/functions/strftime.html
      // TODO: Implement.
      return 0;
    }var _strftime_l=_strftime;
  function __parseInt64(str, endptr, base, min, max, unsign) {
      var isNegative = false;
      // Skip space.
      while (_isspace(HEAP8[(str)])) str++;
      // Check for a plus/minus sign.
      if (HEAP8[(str)] == 45) {
        str++;
        isNegative = true;
      } else if (HEAP8[(str)] == 43) {
        str++;
      }
      // Find base.
      var ok = false;
      var finalBase = base;
      if (!finalBase) {
        if (HEAP8[(str)] == 48) {
          if (HEAP8[((str+1)|0)] == 120 ||
              HEAP8[((str+1)|0)] == 88) {
            finalBase = 16;
            str += 2;
          } else {
            finalBase = 8;
            ok = true; // we saw an initial zero, perhaps the entire thing is just "0"
          }
        }
      } else if (finalBase==16) {
        if (HEAP8[(str)] == 48) {
          if (HEAP8[((str+1)|0)] == 120 ||
              HEAP8[((str+1)|0)] == 88) {
            str += 2;
          }
        }
      }
      if (!finalBase) finalBase = 10;
      start = str;
      // Get digits.
      var chr;
      while ((chr = HEAP8[(str)]) != 0) {
        var digit = parseInt(String.fromCharCode(chr), finalBase);
        if (isNaN(digit)) {
          break;
        } else {
          str++;
          ok = true;
        }
      }
      if (!ok) {
        ___setErrNo(ERRNO_CODES.EINVAL);
        return tempRet0 = 0,0;
      }
      // Set end pointer.
      if (endptr) {
        HEAP32[((endptr)>>2)]=str
      }
      try {
        var numberString = isNegative ? '-'+Pointer_stringify(start, str - start) : Pointer_stringify(start, str - start);
        i64Math.fromString(numberString, finalBase, min, max, unsign);
      } catch(e) {
        ___setErrNo(ERRNO_CODES.ERANGE); // not quite correct
      }
      return tempRet0 = HEAP32[(((tempDoublePtr)+(4))>>2)],HEAP32[((tempDoublePtr)>>2)];
    }function _strtoull(str, endptr, base) {
      return __parseInt64(str, endptr, base, 0, '18446744073709551615', true);  // ULONG_MAX.
    }var _strtoull_l=_strtoull;
  function _strtoll(str, endptr, base) {
      return __parseInt64(str, endptr, base, '-9223372036854775808', '9223372036854775807');  // LLONG_MIN, LLONG_MAX.
    }var _strtoll_l=_strtoll;
  function _uselocale(locale) {
      return 0;
    }
  function ___locale_mb_cur_max() { throw '__locale_mb_cur_max not implemented' }
  function _asprintf(s, format, varargs) {
      return _sprintf(-s, format, varargs);
    }function _vasprintf(s, format, va_arg) {
      return _asprintf(s, format, HEAP32[((va_arg)>>2)]);
    }
  function _vsscanf(s, format, va_arg) {
      return _sscanf(s, format, HEAP32[((va_arg)>>2)]);
    }
  var _llvm_memset_p0i8_i64=_memset;
  function _sbrk(bytes) {
      // Implement a Linux-like 'memory area' for our 'process'.
      // Changes the size of the memory area by |bytes|; returns the
      // address of the previous top ('break') of the memory area
      // We control the "dynamic" memory - DYNAMIC_BASE to DYNAMICTOP
      var self = _sbrk;
      if (!self.called) {
        DYNAMICTOP = alignMemoryPage(DYNAMICTOP); // make sure we start out aligned
        self.called = true;
        assert(Runtime.dynamicAlloc);
        self.alloc = Runtime.dynamicAlloc;
        Runtime.dynamicAlloc = function() { abort('cannot dynamically allocate, sbrk now has control') };
      }
      var ret = DYNAMICTOP;
      if (bytes != 0) self.alloc(bytes);
      return ret;  // Previous break location.
    }
  function ___cxa_call_unexpected(exception) {
      Module.printErr('Unexpected exception thrown, this is not properly supported - aborting');
      ABORT = true;
      throw exception;
    }
  var _sqrtf=Math.sqrt;
  function _llvm_lifetime_start() {}
  function _llvm_lifetime_end() {}
  var _fabsf=Math.abs;
  function _fputc(c, stream) {
      // int fputc(int c, FILE *stream);
      // http://pubs.opengroup.org/onlinepubs/000095399/functions/fputc.html
      var chr = unSign(c & 0xFF);
      HEAP8[((_fputc.ret)|0)]=chr
      var ret = _write(stream, _fputc.ret, 1);
      if (ret == -1) {
        if (FS.streams[stream]) FS.streams[stream].error = true;
        return -1;
      } else {
        return chr;
      }
    }function _puts(s) {
      // int puts(const char *s);
      // http://pubs.opengroup.org/onlinepubs/000095399/functions/puts.html
      // NOTE: puts() always writes an extra newline.
      var stdout = HEAP32[((_stdout)>>2)];
      var ret = _fputs(s, stdout);
      if (ret < 0) {
        return ret;
      } else {
        var newlineRet = _fputc(10, stdout);
        return (newlineRet < 0) ? -1 : ret + 1;
      }
    }
  var Browser={mainLoop:{scheduler:null,shouldPause:false,paused:false,queue:[],pause:function () {
          Browser.mainLoop.shouldPause = true;
        },resume:function () {
          if (Browser.mainLoop.paused) {
            Browser.mainLoop.paused = false;
            Browser.mainLoop.scheduler();
          }
          Browser.mainLoop.shouldPause = false;
        },updateStatus:function () {
          if (Module['setStatus']) {
            var message = Module['statusMessage'] || 'Please wait...';
            var remaining = Browser.mainLoop.remainingBlockers;
            var expected = Browser.mainLoop.expectedBlockers;
            if (remaining) {
              if (remaining < expected) {
                Module['setStatus'](message + ' (' + (expected - remaining) + '/' + expected + ')');
              } else {
                Module['setStatus'](message);
              }
            } else {
              Module['setStatus']('');
            }
          }
        }},isFullScreen:false,pointerLock:false,moduleContextCreatedCallbacks:[],workers:[],init:function () {
        if (!Module["preloadPlugins"]) Module["preloadPlugins"] = []; // needs to exist even in workers
        if (Browser.initted || ENVIRONMENT_IS_WORKER) return;
        Browser.initted = true;
        try {
          new Blob();
          Browser.hasBlobConstructor = true;
        } catch(e) {
          Browser.hasBlobConstructor = false;
          console.log("warning: no blob constructor, cannot create blobs with mimetypes");
        }
        Browser.BlobBuilder = typeof MozBlobBuilder != "undefined" ? MozBlobBuilder : (typeof WebKitBlobBuilder != "undefined" ? WebKitBlobBuilder : (!Browser.hasBlobConstructor ? console.log("warning: no BlobBuilder") : null));
        Browser.URLObject = typeof window != "undefined" ? (window.URL ? window.URL : window.webkitURL) : console.log("warning: cannot create object URLs");
        // Support for plugins that can process preloaded files. You can add more of these to
        // your app by creating and appending to Module.preloadPlugins.
        //
        // Each plugin is asked if it can handle a file based on the file's name. If it can,
        // it is given the file's raw data. When it is done, it calls a callback with the file's
        // (possibly modified) data. For example, a plugin might decompress a file, or it
        // might create some side data structure for use later (like an Image element, etc.).
        function getMimetype(name) {
          return {
            'jpg': 'image/jpeg',
            'jpeg': 'image/jpeg',
            'png': 'image/png',
            'bmp': 'image/bmp',
            'ogg': 'audio/ogg',
            'wav': 'audio/wav',
            'mp3': 'audio/mpeg'
          }[name.substr(name.lastIndexOf('.')+1)];
        }
        var imagePlugin = {};
        imagePlugin['canHandle'] = function(name) {
          return !Module.noImageDecoding && /\.(jpg|jpeg|png|bmp)$/.exec(name);
        };
        imagePlugin['handle'] = function(byteArray, name, onload, onerror) {
          var b = null;
          if (Browser.hasBlobConstructor) {
            try {
              b = new Blob([byteArray], { type: getMimetype(name) });
            } catch(e) {
              Runtime.warnOnce('Blob constructor present but fails: ' + e + '; falling back to blob builder');
            }
          }
          if (!b) {
            var bb = new Browser.BlobBuilder();
            bb.append((new Uint8Array(byteArray)).buffer); // we need to pass a buffer, and must copy the array to get the right data range
            b = bb.getBlob();
          }
          var url = Browser.URLObject.createObjectURL(b);
          var img = new Image();
          img.onload = function() {
            assert(img.complete, 'Image ' + name + ' could not be decoded');
            var canvas = document.createElement('canvas');
            canvas.width = img.width;
            canvas.height = img.height;
            var ctx = canvas.getContext('2d');
            ctx.drawImage(img, 0, 0);
            Module["preloadedImages"][name] = canvas;
            Browser.URLObject.revokeObjectURL(url);
            if (onload) onload(byteArray);
          };
          img.onerror = function(event) {
            console.log('Image ' + url + ' could not be decoded');
            if (onerror) onerror();
          };
          img.src = url;
        };
        Module['preloadPlugins'].push(imagePlugin);
        var audioPlugin = {};
        audioPlugin['canHandle'] = function(name) {
          return !Module.noAudioDecoding && name.substr(-4) in { '.ogg': 1, '.wav': 1, '.mp3': 1 };
        };
        audioPlugin['handle'] = function(byteArray, name, onload, onerror) {
          var done = false;
          function finish(audio) {
            if (done) return;
            done = true;
            Module["preloadedAudios"][name] = audio;
            if (onload) onload(byteArray);
          }
          function fail() {
            if (done) return;
            done = true;
            Module["preloadedAudios"][name] = new Audio(); // empty shim
            if (onerror) onerror();
          }
          if (Browser.hasBlobConstructor) {
            try {
              var b = new Blob([byteArray], { type: getMimetype(name) });
            } catch(e) {
              return fail();
            }
            var url = Browser.URLObject.createObjectURL(b); // XXX we never revoke this!
            var audio = new Audio();
            audio.addEventListener('canplaythrough', function() { finish(audio) }, false); // use addEventListener due to chromium bug 124926
            audio.onerror = function(event) {
              if (done) return;
              console.log('warning: browser could not fully decode audio ' + name + ', trying slower base64 approach');
              function encode64(data) {
                var BASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/';
                var PAD = '=';
                var ret = '';
                var leftchar = 0;
                var leftbits = 0;
                for (var i = 0; i < data.length; i++) {
                  leftchar = (leftchar << 8) | data[i];
                  leftbits += 8;
                  while (leftbits >= 6) {
                    var curr = (leftchar >> (leftbits-6)) & 0x3f;
                    leftbits -= 6;
                    ret += BASE[curr];
                  }
                }
                if (leftbits == 2) {
                  ret += BASE[(leftchar&3) << 4];
                  ret += PAD + PAD;
                } else if (leftbits == 4) {
                  ret += BASE[(leftchar&0xf) << 2];
                  ret += PAD;
                }
                return ret;
              }
              audio.src = 'data:audio/x-' + name.substr(-3) + ';base64,' + encode64(byteArray);
              finish(audio); // we don't wait for confirmation this worked - but it's worth trying
            };
            audio.src = url;
            // workaround for chrome bug 124926 - we do not always get oncanplaythrough or onerror
            Browser.safeSetTimeout(function() {
              finish(audio); // try to use it even though it is not necessarily ready to play
            }, 10000);
          } else {
            return fail();
          }
        };
        Module['preloadPlugins'].push(audioPlugin);
        // Canvas event setup
        var canvas = Module['canvas'];
        canvas.requestPointerLock = canvas['requestPointerLock'] ||
                                    canvas['mozRequestPointerLock'] ||
                                    canvas['webkitRequestPointerLock'];
        canvas.exitPointerLock = document['exitPointerLock'] ||
                                 document['mozExitPointerLock'] ||
                                 document['webkitExitPointerLock'] ||
                                 function(){}; // no-op if function does not exist
        canvas.exitPointerLock = canvas.exitPointerLock.bind(document);
        function pointerLockChange() {
          Browser.pointerLock = document['pointerLockElement'] === canvas ||
                                document['mozPointerLockElement'] === canvas ||
                                document['webkitPointerLockElement'] === canvas;
        }
        document.addEventListener('pointerlockchange', pointerLockChange, false);
        document.addEventListener('mozpointerlockchange', pointerLockChange, false);
        document.addEventListener('webkitpointerlockchange', pointerLockChange, false);
        if (Module['elementPointerLock']) {
          canvas.addEventListener("click", function(ev) {
            if (!Browser.pointerLock && canvas.requestPointerLock) {
              canvas.requestPointerLock();
              ev.preventDefault();
            }
          }, false);
        }
      },createContext:function (canvas, useWebGL, setInModule) {
        var ctx;
        try {
          if (useWebGL) {
            ctx = canvas.getContext('experimental-webgl', {
              alpha: false
            });
          } else {
            ctx = canvas.getContext('2d');
          }
          if (!ctx) throw ':(';
        } catch (e) {
          Module.print('Could not create canvas - ' + e);
          return null;
        }
        if (useWebGL) {
          // Set the background of the WebGL canvas to black
          canvas.style.backgroundColor = "black";
          // Warn on context loss
          canvas.addEventListener('webglcontextlost', function(event) {
            alert('WebGL context lost. You will need to reload the page.');
          }, false);
        }
        if (setInModule) {
          Module.ctx = ctx;
          Module.useWebGL = useWebGL;
          Browser.moduleContextCreatedCallbacks.forEach(function(callback) { callback() });
          Browser.init();
        }
        return ctx;
      },destroyContext:function (canvas, useWebGL, setInModule) {},fullScreenHandlersInstalled:false,lockPointer:undefined,resizeCanvas:undefined,requestFullScreen:function (lockPointer, resizeCanvas) {
        Browser.lockPointer = lockPointer;
        Browser.resizeCanvas = resizeCanvas;
        if (typeof Browser.lockPointer === 'undefined') Browser.lockPointer = true;
        if (typeof Browser.resizeCanvas === 'undefined') Browser.resizeCanvas = false;
        var canvas = Module['canvas'];
        function fullScreenChange() {
          Browser.isFullScreen = false;
          if ((document['webkitFullScreenElement'] || document['webkitFullscreenElement'] ||
               document['mozFullScreenElement'] || document['mozFullscreenElement'] ||
               document['fullScreenElement'] || document['fullscreenElement']) === canvas) {
            canvas.cancelFullScreen = document['cancelFullScreen'] ||
                                      document['mozCancelFullScreen'] ||
                                      document['webkitCancelFullScreen'];
            canvas.cancelFullScreen = canvas.cancelFullScreen.bind(document);
            if (Browser.lockPointer) canvas.requestPointerLock();
            Browser.isFullScreen = true;
            if (Browser.resizeCanvas) Browser.setFullScreenCanvasSize();
          } else if (Browser.resizeCanvas){
            Browser.setWindowedCanvasSize();
          }
          if (Module['onFullScreen']) Module['onFullScreen'](Browser.isFullScreen);
        }
        if (!Browser.fullScreenHandlersInstalled) {
          Browser.fullScreenHandlersInstalled = true;
          document.addEventListener('fullscreenchange', fullScreenChange, false);
          document.addEventListener('mozfullscreenchange', fullScreenChange, false);
          document.addEventListener('webkitfullscreenchange', fullScreenChange, false);
        }
        canvas.requestFullScreen = canvas['requestFullScreen'] ||
                                   canvas['mozRequestFullScreen'] ||
                                   (canvas['webkitRequestFullScreen'] ? function() { canvas['webkitRequestFullScreen'](Element['ALLOW_KEYBOARD_INPUT']) } : null);
        canvas.requestFullScreen();
      },requestAnimationFrame:function (func) {
        if (!window.requestAnimationFrame) {
          window.requestAnimationFrame = window['requestAnimationFrame'] ||
                                         window['mozRequestAnimationFrame'] ||
                                         window['webkitRequestAnimationFrame'] ||
                                         window['msRequestAnimationFrame'] ||
                                         window['oRequestAnimationFrame'] ||
                                         window['setTimeout'];
        }
        window.requestAnimationFrame(func);
      },safeCallback:function (func) {
        return function() {
          if (!ABORT) return func.apply(null, arguments);
        };
      },safeRequestAnimationFrame:function (func) {
        return Browser.requestAnimationFrame(function() {
          if (!ABORT) func();
        });
      },safeSetTimeout:function (func, timeout) {
        return setTimeout(function() {
          if (!ABORT) func();
        }, timeout);
      },safeSetInterval:function (func, timeout) {
        return setInterval(function() {
          if (!ABORT) func();
        }, timeout);
      },getUserMedia:function (func) {
        if(!window.getUserMedia) {
          window.getUserMedia = navigator['getUserMedia'] ||
                                navigator['mozGetUserMedia'];
        }
        window.getUserMedia(func);
      },getMovementX:function (event) {
        return event['movementX'] ||
               event['mozMovementX'] ||
               event['webkitMovementX'] ||
               0;
      },getMovementY:function (event) {
        return event['movementY'] ||
               event['mozMovementY'] ||
               event['webkitMovementY'] ||
               0;
      },mouseX:0,mouseY:0,mouseMovementX:0,mouseMovementY:0,calculateMouseEvent:function (event) { // event should be mousemove, mousedown or mouseup
        if (Browser.pointerLock) {
          // When the pointer is locked, calculate the coordinates
          // based on the movement of the mouse.
          // Workaround for Firefox bug 764498
          if (event.type != 'mousemove' &&
              ('mozMovementX' in event)) {
            Browser.mouseMovementX = Browser.mouseMovementY = 0;
          } else {
            Browser.mouseMovementX = Browser.getMovementX(event);
            Browser.mouseMovementY = Browser.getMovementY(event);
          }
          Browser.mouseX = SDL.mouseX + Browser.mouseMovementX;
          Browser.mouseY = SDL.mouseY + Browser.mouseMovementY;
        } else {
          // Otherwise, calculate the movement based on the changes
          // in the coordinates.
          var rect = Module["canvas"].getBoundingClientRect();
          var x = event.pageX - (window.scrollX + rect.left);
          var y = event.pageY - (window.scrollY + rect.top);
          // the canvas might be CSS-scaled compared to its backbuffer;
          // SDL-using content will want mouse coordinates in terms
          // of backbuffer units.
          var cw = Module["canvas"].width;
          var ch = Module["canvas"].height;
          x = x * (cw / rect.width);
          y = y * (ch / rect.height);
          Browser.mouseMovementX = x - Browser.mouseX;
          Browser.mouseMovementY = y - Browser.mouseY;
          Browser.mouseX = x;
          Browser.mouseY = y;
        }
      },xhrLoad:function (url, onload, onerror) {
        var xhr = new XMLHttpRequest();
        xhr.open('GET', url, true);
        xhr.responseType = 'arraybuffer';
        xhr.onload = function() {
          if (xhr.status == 200 || (xhr.status == 0 && xhr.response)) { // file URLs can return 0
            onload(xhr.response);
          } else {
            onerror();
          }
        };
        xhr.onerror = onerror;
        xhr.send(null);
      },asyncLoad:function (url, onload, onerror, noRunDep) {
        Browser.xhrLoad(url, function(arrayBuffer) {
          assert(arrayBuffer, 'Loading data file "' + url + '" failed (no arrayBuffer).');
          onload(new Uint8Array(arrayBuffer));
          if (!noRunDep) removeRunDependency('al ' + url);
        }, function(event) {
          if (onerror) {
            onerror();
          } else {
            throw 'Loading data file "' + url + '" failed.';
          }
        });
        if (!noRunDep) addRunDependency('al ' + url);
      },resizeListeners:[],updateResizeListeners:function () {
        var canvas = Module['canvas'];
        Browser.resizeListeners.forEach(function(listener) {
          listener(canvas.width, canvas.height);
        });
      },setCanvasSize:function (width, height, noUpdates) {
        var canvas = Module['canvas'];
        canvas.width = width;
        canvas.height = height;
        if (!noUpdates) Browser.updateResizeListeners();
      },windowedWidth:0,windowedHeight:0,setFullScreenCanvasSize:function () {
        var canvas = Module['canvas'];
        this.windowedWidth = canvas.width;
        this.windowedHeight = canvas.height;
        canvas.width = screen.width;
        canvas.height = screen.height;
        var flags = HEAPU32[((SDL.screen+Runtime.QUANTUM_SIZE*0)>>2)];
        flags = flags | 0x00800000; // set SDL_FULLSCREEN flag
        HEAP32[((SDL.screen+Runtime.QUANTUM_SIZE*0)>>2)]=flags
        Browser.updateResizeListeners();
      },setWindowedCanvasSize:function () {
        var canvas = Module['canvas'];
        canvas.width = this.windowedWidth;
        canvas.height = this.windowedHeight;
        var flags = HEAPU32[((SDL.screen+Runtime.QUANTUM_SIZE*0)>>2)];
        flags = flags & ~0x00800000; // clear SDL_FULLSCREEN flag
        HEAP32[((SDL.screen+Runtime.QUANTUM_SIZE*0)>>2)]=flags
        Browser.updateResizeListeners();
      }};
__ATINIT__.unshift({ func: function() { if (!Module["noFSInit"] && !FS.init.initialized) FS.init() } });__ATMAIN__.push({ func: function() { FS.ignorePermissions = false } });__ATEXIT__.push({ func: function() { FS.quit() } });Module["FS_createFolder"] = FS.createFolder;Module["FS_createPath"] = FS.createPath;Module["FS_createDataFile"] = FS.createDataFile;Module["FS_createPreloadedFile"] = FS.createPreloadedFile;Module["FS_createLazyFile"] = FS.createLazyFile;Module["FS_createLink"] = FS.createLink;Module["FS_createDevice"] = FS.createDevice;
___errno_state = Runtime.staticAlloc(4); HEAP32[((___errno_state)>>2)]=0;
_fgetc.ret = allocate([0], "i8", ALLOC_STATIC);
___buildEnvironment(ENV);
___strtok_state = Runtime.staticAlloc(4);
_llvm_eh_exception.buf = allocate(12, "void*", ALLOC_STATIC);
_fputc.ret = allocate([0], "i8", ALLOC_STATIC);
Module["requestFullScreen"] = function(lockPointer, resizeCanvas) { Browser.requestFullScreen(lockPointer, resizeCanvas) };
  Module["requestAnimationFrame"] = function(func) { Browser.requestAnimationFrame(func) };
  Module["pauseMainLoop"] = function() { Browser.mainLoop.pause() };
  Module["resumeMainLoop"] = function() { Browser.mainLoop.resume() };
  Module["getUserMedia"] = function() { Browser.getUserMedia() }
STACK_BASE = STACKTOP = Runtime.alignMemory(STATICTOP);
staticSealed = true; // seal the static portion of memory
STACK_MAX = STACK_BASE + 5242880;
DYNAMIC_BASE = DYNAMICTOP = Runtime.alignMemory(STACK_MAX);
assert(DYNAMIC_BASE < TOTAL_MEMORY); // Stack must fit in TOTAL_MEMORY; allocations from here on may enlarge TOTAL_MEMORY
var FUNCTION_TABLE = [0,0,__ZNSt3__18messagesIwED0Ev,0,_subgraph_search,0,_spline_merge4418,0,_arrow_type_box,0,__ZNSt3__17num_getIcNS_19istreambuf_iteratorIcNS_11char_traitsIcEEEEED0Ev
,0,__ZNKSt3__18numpunctIcE12do_falsenameEv,0,__Z9agedgegetP8Agedge_sRKNSt3__112basic_stringIcNS1_11char_traitsIcEENS1_9allocatorIcEEEE,0,_makeAddPoly,0,_xdot_polygon,0,___ZN10emscripten8internal7InvokerINSt3__112basic_stringIcNS2_11char_traitsIcEENS2_9allocatorIcEEEEJP8Agnode_sRKS8_EE6invokeEPFS8_SA_SC_ESA_PNS0_11BindingTypeIS8_E3$_0E_
,0,__ZN10emscripten8internal7InvokerIP8Agedge_sJP8Agraph_sP8Agnode_sEE6invokeEPFS3_S5_S7_ES5_S7_,0,_agdictobjfree,0,__ZNKSt3__120__time_get_c_storageIwE3__rEv,0,__ZNSt3__18messagesIcED0Ev,0,_epsf_gencode
,0,_zoom_in_cb,0,__ZN10emscripten8internal13getActualTypeI8Agedge_sEEPKNS0_7_TYPEIDEPT_,0,__ZNKSt3__18time_getIcNS_19istreambuf_iteratorIcNS_11char_traitsIcEEEEE11do_get_timeES4_S4_RNS_8ios_baseERjP2tm,0,_sfdp_cleanup,0,__ZNKSt3__18time_getIcNS_19istreambuf_iteratorIcNS_11char_traitsIcEEEEE11do_get_yearES4_S4_RNS_8ios_baseERjP2tm
,0,_psgen_begin_graph,0,__ZNSt12length_errorD0Ev,0,_arrow_type_dot,0,_sfdp_layout,0,_point_init
,0,__ZNSt3__115basic_streambufIcNS_11char_traitsIcEEED1Ev,0,_pov_begin_node,0,__ZNKSt3__18time_getIwNS_19istreambuf_iteratorIwNS_11char_traitsIwEEEEE11do_get_timeES4_S4_RNS_8ios_baseERjP2tm,0,_xdot_end_cluster,0,_psgen_bezier
,0,__ZNKSt3__15ctypeIcE10do_toupperEc,0,__ZNSt3__16locale2id6__initEv,0,_ordercmpf,0,__ZNKSt3__18time_getIwNS_19istreambuf_iteratorIwNS_11char_traitsIwEEEEE14do_get_weekdayES4_S4_RNS_8ios_baseERjP2tm,0,__ZNKSt3__17codecvtIDsc10_mbstate_tE9do_lengthERS1_PKcS5_j
,0,__ZNSt3__110__stdinbufIcE9pbackfailEi,0,__ZNSt3__115basic_streambufIcNS_11char_traitsIcEEE9underflowEv,0,__ZNSt3__110__stdinbufIwED0Ev,0,_svg_begin_job,0,__ZNSt11logic_errorD2Ev
,0,_fig_begin_edge,0,___ZN10emscripten8internal7InvokerINSt3__112basic_stringIcNS2_11char_traitsIcEENS2_9allocatorIcEEEEJP5GVC_sP8Agraph_sRKS8_EE6invokeEPFS8_SA_SC_SE_ESA_SC_PNS0_11BindingTypeIS8_E3$_0E_,0,__ZNKSt3__17collateIcE7do_hashEPKcS3_,0,_cmpDegree,0,_subedge_search
,0,__ZNKSt3__120__time_get_c_storageIwE8__monthsEv,0,__ZNKSt3__18time_getIcNS_19istreambuf_iteratorIcNS_11char_traitsIcEEEEE6do_getES4_S4_RNS_8ios_baseERjP2tmcc,0,_pov_begin_cluster,0,__ZNSt3__17num_getIwNS_19istreambuf_iteratorIwNS_11char_traitsIwEEEEED0Ev,0,__ZNKSt3__19money_putIwNS_19ostreambuf_iteratorIwNS_11char_traitsIwEEEEE6do_putES4_bRNS_8ios_baseEwRKNS_12basic_stringIwS3_NS_9allocatorIwEEEE
,0,_psgen_polygon,0,___ZN10emscripten8internal7InvokerIiJP8Agraph_sRKNSt3__112basic_stringIcNS4_11char_traitsIcEENS4_9allocatorIcEEEESC_EE6invokeEPFiS3_SC_SC_ES3_PNS0_11BindingTypeISA_E3$_0ESJ__,0,__Z10aggraphsetP8Agraph_sRKNSt3__112basic_stringIcNS1_11char_traitsIcEENS1_9allocatorIcEEEES9_,0,_cmppair,0,__ZNSt12out_of_rangeD0Ev
,0,__ZNKSt3__17codecvtIcc10_mbstate_tE6do_outERS1_PKcS5_RS5_PcS7_RS7_,0,___ZN10emscripten8internal7InvokerINSt3__112basic_stringIcNS2_11char_traitsIcEENS2_9allocatorIcEEEEJP8Agedge_sRKS8_EE6invokeEPFS8_SA_SC_ESA_PNS0_11BindingTypeIS8_E3$_0E_,0,__ZNSt3__19money_getIwNS_19istreambuf_iteratorIwNS_11char_traitsIwEEEEED0Ev,0,_ecmp,0,_pov_end_cluster
,0,__ZNKSt3__110moneypunctIwLb1EE16do_positive_signEv,0,_vml_ellipse,0,__ZNKSt3__15ctypeIwE10do_tolowerEPwPKw,0,__ZNSt3__115basic_streambufIcNS_11char_traitsIcEEE5uflowEv,0,_gvevent_button_press
,0,_my_init_graph,0,_psgen_ellipse,0,__ZNSt3__17collateIcED1Ev,0,__ZNSt3__18ios_base7failureD2Ev,0,_pov_end_edge
,0,_psgen_end_page,0,_tkgen_polyline,0,__ZNK10__cxxabiv121__vmi_class_type_info16search_above_dstEPNS_19__dynamic_cast_infoEPKvS4_ib,0,_arrow_type_diamond,0,__ZNSt9bad_allocD2Ev
,0,_fcmpf,0,_svg_end_page,0,_gvrender_comparestr,0,__Z13wrap_gvLayoutP5GVC_sP8Agraph_sRKNSt3__112basic_stringIcNS3_11char_traitsIcEENS3_9allocatorIcEEEE,0,_psgen_end_edge
,0,__ZNKSt3__17codecvtIDsc10_mbstate_tE10do_unshiftERS1_PcS4_RS4_,0,_pic_end_graph,0,__ZNSt3__16locale5facetD0Ev,0,_psgen_textpara,0,_psgen_begin_node
,0,__ZNKSt3__120__time_get_c_storageIwE3__cEv,0,__ZNKSt3__17num_putIwNS_19ostreambuf_iteratorIwNS_11char_traitsIwEEEEE6do_putES4_RNS_8ios_baseEwy,0,__ZNKSt3__17num_putIwNS_19ostreambuf_iteratorIwNS_11char_traitsIwEEEEE6do_putES4_RNS_8ios_baseEwx,0,__ZNSt3__15ctypeIcED0Ev,0,_memresize
,0,_arrow_type_curve,0,__ZNKSt3__17num_putIwNS_19ostreambuf_iteratorIwNS_11char_traitsIwEEEEE6do_putES4_RNS_8ios_baseEwm,0,__ZNSt3__111__stdoutbufIcE4syncEv,0,_psgen_end_cluster,0,__ZNSt3__18time_getIwNS_19istreambuf_iteratorIwNS_11char_traitsIwEEEEED1Ev
,0,__ZNSt3__17num_getIwNS_19istreambuf_iteratorIwNS_11char_traitsIwEEEEED1Ev,0,__ZNKSt3__17num_putIwNS_19ostreambuf_iteratorIwNS_11char_traitsIwEEEEE6do_putES4_RNS_8ios_baseEwe,0,__ZNKSt3__17num_putIwNS_19ostreambuf_iteratorIwNS_11char_traitsIwEEEEE6do_putES4_RNS_8ios_baseEwd,0,__ZNKSt3__110moneypunctIcLb1EE16do_decimal_pointEv,0,__ZNKSt3__17codecvtIwc10_mbstate_tE11do_encodingEv
,0,__ZNKSt3__17num_putIwNS_19ostreambuf_iteratorIwNS_11char_traitsIwEEEEE6do_putES4_RNS_8ios_baseEwb,0,_freef,0,_simple_delrec,0,_agedgeseqcmpf,0,__ZNK10__cxxabiv120__si_class_type_info16search_above_dstEPNS_19__dynamic_cast_infoEPKvS4_ib
,0,__ZNKSt3__19money_putIcNS_19ostreambuf_iteratorIcNS_11char_traitsIcEEEEE6do_putES4_bRNS_8ios_baseEcRKNS_12basic_stringIcS3_NS_9allocatorIcEEEE,0,_dot_cleanup,0,_map_begin_page,0,__Z11wrap_agreadRKNSt3__112basic_stringIcNS_11char_traitsIcEENS_9allocatorIcEEEE,0,__ZNSt3__113basic_ostreamIwNS_11char_traitsIwEEED1Ev
,0,_fig_comment,0,_agdeledge,0,__ZNKSt3__17codecvtIwc10_mbstate_tE9do_lengthERS1_PKcS5_j,0,_gvputs,0,__ZNSt3__19money_putIcNS_19ostreambuf_iteratorIcNS_11char_traitsIcEEEEED1Ev
,0,__ZNKSt3__18time_getIcNS_19istreambuf_iteratorIcNS_11char_traitsIcEEEEE13do_date_orderEv,0,__ZNSt3__18messagesIcED1Ev,0,_freeItem3766,0,__ZNKSt3__120__time_get_c_storageIwE7__weeksEv,0,_osage_cleanup
,0,_fig_end_graph,0,_record_free,0,__ZNSt3__16locale5facet16__on_zero_sharedEv,0,__ZN10emscripten8internal7InvokerIP8Agnode_sJP8Agraph_smiEE6invokeEPFS3_S5_miES5_mi,0,_arrow_type_normal
,0,__ZNKSt3__15ctypeIwE8do_widenEc,0,__ZNKSt3__18time_putIwNS_19ostreambuf_iteratorIwNS_11char_traitsIwEEEEE6do_putES4_RNS_8ios_baseEwPK2tmcc,0,_ps_image_free,0,__ZNSt3__110__stdinbufIcE5uflowEv,0,_pov_end_page
,0,_gt,0,__ZNSt3__17num_putIwNS_19ostreambuf_iteratorIwNS_11char_traitsIwEEEEED1Ev,0,_svg_begin_cluster,0,__ZNSt3__119__iostream_categoryD0Ev,0,__ZTv0_n12_NSt3__113basic_istreamIcNS_11char_traitsIcEEED0Ev
,0,__ZNSt3__17num_getIcNS_19istreambuf_iteratorIcNS_11char_traitsIcEEEEED1Ev,0,__ZNSt3__115basic_streambufIwNS_11char_traitsIwEEE5uflowEv,0,_memfree,0,__ZNKSt3__110moneypunctIwLb0EE13do_neg_formatEv,0,_cmpItem
,0,_objdelrec,0,_xdot_end_node,0,_svg_end_cluster,0,_dot_begin_graph,0,__ZNKSt3__17codecvtIwc10_mbstate_tE5do_inERS1_PKcS5_RS5_PwS7_RS7_
,0,__ZN10emscripten8internal7InvokerINSt3__112basic_stringIcNS2_11char_traitsIcEENS2_9allocatorIcEEEEJP5GVC_sEE6invokeEPFS8_SA_ESA_,0,__ZNKSt3__17codecvtIDsc10_mbstate_tE5do_inERS1_PKcS5_RS5_PDsS7_RS7_,0,__ZNKSt3__15ctypeIcE8do_widenEc,0,__ZNSt3__110moneypunctIwLb0EED0Ev,0,__ZNKSt3__18time_getIwNS_19istreambuf_iteratorIwNS_11char_traitsIwEEEEE13do_date_orderEv
,0,__ZNKSt3__17codecvtIDic10_mbstate_tE9do_lengthERS1_PKcS5_j,0,_svg_begin_layer,0,_vml_textpara,0,__ZNSt3__16locale5__impD2Ev,0,_fig_polyline
,0,__ZNSt3__115basic_streambufIwNS_11char_traitsIwEEE9underflowEv,0,_agsubnodeidcmpf,0,_zoom_out_cb,0,_pic_begin_page,0,_patchwork_cleanup
,0,_colorcmpf,0,_record_port,0,_xdot_textpara,0,_freeitem,0,__ZNKSt3__17num_putIwNS_19ostreambuf_iteratorIwNS_11char_traitsIwEEEEE6do_putES4_RNS_8ios_baseEwPKv
,0,_circo_layout,0,__ZNSt3__18numpunctIcED2Ev,0,_vml_begin_job,0,_osage_layout,0,__ZNSt3__115basic_streambufIcNS_11char_traitsIcEEE8overflowEi
,0,_xdot_bezier,0,__ZNSt3__17codecvtIcc10_mbstate_tED0Ev,0,_agdeledgeimage,0,__ZNKSt3__18numpunctIcE11do_groupingEv,0,_map_begin_anchor
,0,__ZNSt3__115basic_streambufIwNS_11char_traitsIwEEE4syncEv,0,__ZN10__cxxabiv119__pointer_type_infoD0Ev,0,_epsf_init,0,_pic_textpara,0,__ZNSt3__110__stdinbufIcED1Ev
,0,_point_inside,0,__ZNKSt3__120__time_get_c_storageIwE3__xEv,0,__ZNKSt3__110moneypunctIcLb1EE13do_neg_formatEv,0,_dtvsearch,0,__ZNSt3__110__stdinbufIwE9pbackfailEi
,0,_nodeposcmpf,0,_gvevent_layout,0,__ZNSt3__18numpunctIcED0Ev,0,__ZNSt3__111__stdoutbufIcE8overflowEi,0,_pic_end_page
,0,_twopi_cleanup,0,__ZNSt3__119__iostream_categoryD1Ev,0,__ZNKSt3__120__time_get_c_storageIwE7__am_pmEv,0,__ZTv0_n12_NSt3__113basic_ostreamIwNS_11char_traitsIwEEED1Ev,0,__ZN10emscripten8internal14raw_destructorI8Agraph_sEEvPT_
,0,__ZNSt3__111__stdoutbufIwE5imbueERKNS_6localeE,0,_iofread,0,_sortf,0,__ZNKSt3__18messagesIcE8do_closeEi,0,_fig_bezier
,0,__ZNKSt3__15ctypeIwE5do_isEPKwS3_Pt,0,__ZNSt13runtime_errorD2Ev,0,__ZNKSt3__15ctypeIwE10do_toupperEw,0,__ZNKSt3__110moneypunctIwLb1EE16do_negative_signEv,0,__ZNKSt3__15ctypeIwE9do_narrowEPKwS3_cPc
,0,__ZNKSt3__17codecvtIDic10_mbstate_tE11do_encodingEv,0,_fdp_layout,0,__Z14wrap_gvContextv,0,__ZNSt3__115basic_streambufIwNS_11char_traitsIwEEE5imbueERKNS_6localeE,0,__ZNKSt3__110moneypunctIcLb0EE16do_negative_signEv
,0,_pic_begin_graph,0,_record_init,0,__ZNSt3__17collateIwED1Ev,0,__ZNKSt3__18time_getIwNS_19istreambuf_iteratorIwNS_11char_traitsIwEEEEE16do_get_monthnameES4_S4_RNS_8ios_baseERjP2tm,0,__ZNK10__cxxabiv117__class_type_info9can_catchEPKNS_16__shim_type_infoERPv
,0,__ZNKSt8bad_cast4whatEv,0,_agsubnodeseqcmpf,0,__ZNSt3__110moneypunctIcLb0EED1Ev,0,__ZNKSt3__18messagesIcE6do_getEiiiRKNS_12basic_stringIcNS_11char_traitsIcEENS_9allocatorIcEEEE,0,_edgelblcmpfn
,0,__ZNSt3__18numpunctIwED2Ev,0,_pov_begin_job,0,_svg_end_layer,0,_my_init_node,0,_tkgen_begin_edge
,0,_pov_begin_graph,0,__ZNKSt13runtime_error4whatEv,0,_vml_end_graph,0,_core_loadimage_vrml,0,_poly_path
,0,__ZNSt3__19money_putIwNS_19ostreambuf_iteratorIwNS_11char_traitsIwEEEEED0Ev,0,_compoundEdges,0,_newitem4397,0,__ZNSt3__117__widen_from_utf8ILj32EED0Ev,0,_core_loadimage_vml
,0,__ZNKSt3__18numpunctIwE16do_thousands_sepEv,0,_down_cb,0,_agprvnode,0,__ZNKSt3__18time_getIwNS_19istreambuf_iteratorIwNS_11char_traitsIwEEEEE6do_getES4_S4_RNS_8ios_baseERjP2tmcc,0,__ZNSt3__113basic_istreamIwNS_11char_traitsIwEEED1Ev
,0,_svg_begin_node,0,__ZNKSt3__18numpunctIcE16do_decimal_pointEv,0,_vml_polyline,0,__ZN10emscripten8internal7InvokerIiJP5GVC_sP8Agraph_sEE6invokeEPFiS3_S5_ES3_S5_,0,_svg_begin_edge
,0,__Z10aggraphgetP8Agraph_sRKNSt3__112basic_stringIcNS1_11char_traitsIcEENS1_9allocatorIcEEEE,0,_tkgen_begin_node,0,__ZNKSt3__110moneypunctIwLb0EE16do_negative_signEv,0,__ZNK10__cxxabiv120__si_class_type_info16search_below_dstEPNS_19__dynamic_cast_infoEPKvib,0,_agnode
,0,__ZNKSt3__120__time_get_c_storageIcE3__xEv,0,__ZNSt3__17collateIwED0Ev,0,__ZNKSt3__110moneypunctIcLb0EE16do_positive_signEv,0,__ZNSt3__18time_putIwNS_19ostreambuf_iteratorIwNS_11char_traitsIwEEEEED0Ev,0,_agfstnode
,0,__ZNKSt3__17codecvtIwc10_mbstate_tE13do_max_lengthEv,0,__ZNSt3__17num_putIcNS_19ostreambuf_iteratorIcNS_11char_traitsIcEEEEED0Ev,0,_dtmemory,0,_svg_begin_page,0,__ZNK10__cxxabiv116__shim_type_info5noop2Ev
,0,_nodecmp,0,___ZN10emscripten8internal7InvokerIP8Agraph_sJRKNSt3__112basic_stringIcNS4_11char_traitsIcEENS4_9allocatorIcEEEEEE6invokeEPFS3_SC_EPNS0_11BindingTypeISA_E3$_0E_,0,__ZNSt3__17num_putIcNS_19ostreambuf_iteratorIcNS_11char_traitsIcEEEEED1Ev,0,__ZNK10__cxxabiv116__shim_type_info5noop1Ev,0,___ZN10emscripten8internal7InvokerIiJP5GVC_sP8Agraph_sRKNSt3__112basic_stringIcNS6_11char_traitsIcEENS6_9allocatorIcEEEEEE6invokeEPFiS3_S5_SE_ES3_S5_PNS0_11BindingTypeISC_E3$_0E_
,0,__ZNKSt3__18numpunctIwE16do_decimal_pointEv,0,_pov_polygon,0,__ZN10emscripten8internal13getActualTypeI8Agraph_sEEPKNS0_7_TYPEIDEPT_,0,__ZNSt3__115basic_streambufIcNS_11char_traitsIcEEE4syncEv,0,_psgen_library_shape
,0,__ZNK10__cxxabiv117__class_type_info16search_below_dstEPNS_19__dynamic_cast_infoEPKvib,0,_fig_end_node,0,__ZNKSt3__110moneypunctIcLb0EE11do_groupingEv,0,__ZNK10__cxxabiv120__si_class_type_info27has_unambiguous_public_baseEPNS_19__dynamic_cast_infoEPvi,0,_fig_begin_page
,0,__ZNKSt3__110moneypunctIwLb1EE14do_frac_digitsEv,0,_agnedges,0,__ZNK10__cxxabiv121__vmi_class_type_info27has_unambiguous_public_baseEPNS_19__dynamic_cast_infoEPvi,0,__ZNKSt3__17num_putIwNS_19ostreambuf_iteratorIwNS_11char_traitsIwEEEEE6do_putES4_RNS_8ios_baseEwl,0,__ZNKSt3__120__time_get_c_storageIcE3__XEv
,0,__ZNKSt3__17num_putIcNS_19ostreambuf_iteratorIcNS_11char_traitsIcEEEEE6do_putES4_RNS_8ios_baseEcm,0,__ZNKSt3__15ctypeIwE9do_narrowEwc,0,__ZN10__cxxabiv123__fundamental_type_infoD0Ev,0,__Z9agnodegetP8Agnode_sRKNSt3__112basic_stringIcNS1_11char_traitsIcEENS1_9allocatorIcEEEE,0,_neato_layout
,0,__ZN10emscripten8internal15raw_constructorI8Agdesc_sJEEEPT_DpNS0_11BindingTypeIT0_E8WireTypeE,0,_pov_textpara,0,_cmpitems,0,___ZN10emscripten8internal7InvokerIiJP8Agedge_sRKNSt3__112basic_stringIcNS4_11char_traitsIcEENS4_9allocatorIcEEEESC_EE6invokeEPFiS3_SC_SC_ES3_PNS0_11BindingTypeISA_E3$_0ESJ__,0,__ZNSt3__111__stdoutbufIwE4syncEv
,0,__ZN10emscripten8internal14raw_destructorI5GVC_sEEvPT_,0,_pov_begin_layer,0,__ZNSt3__110moneypunctIwLb0EED1Ev,0,__ZN10emscripten8internal7InvokerIiJP5GVC_sEE6invokeEPFiS3_ES3_,0,_core_loadimage_null
,0,_agnnodes,0,__ZNKSt3__17codecvtIcc10_mbstate_tE10do_unshiftERS1_PcS4_RS4_,0,_right_cb,0,__ZNSt3__113basic_istreamIcNS_11char_traitsIcEEED1Ev,0,_mkItem
,0,_xdot_ellipse,0,_pic_comment,0,__ZTv0_n12_NSt3__113basic_ostreamIcNS_11char_traitsIcEEED1Ev,0,_svg_polyline,0,_poly_init
,0,__ZNSt3__18time_putIcNS_19ostreambuf_iteratorIcNS_11char_traitsIcEEEEED0Ev,0,__ZNSt3__18time_getIwNS_19istreambuf_iteratorIwNS_11char_traitsIwEEEEED0Ev,0,__ZNKSt3__17collateIwE7do_hashEPKwS3_,0,_agedgeidcmpf,0,__ZNSt3__115basic_streambufIcNS_11char_traitsIcEEE7seekposENS_4fposI10_mbstate_tEEj
,0,__ZNSt3__111__stdoutbufIcE5imbueERKNS_6localeE,0,_agdegree,0,_core_loadimage_svg,0,__ZNKSt3__110moneypunctIcLb1EE16do_thousands_sepEv,0,_pov_begin_edge
,0,__ZNSt3__18ios_baseD0Ev,0,__ZNSt3__111__stdoutbufIcED1Ev,0,__ZNSt3__110moneypunctIcLb1EED0Ev,0,_fig_begin_graph,0,__ZNSt9bad_allocD0Ev
,0,__ZNSt3__115basic_streambufIcNS_11char_traitsIcEEED0Ev,0,_agdelnodeimage,0,_xdot_end_edge,0,_agfstin,0,_Operator_diag_precon_apply
,0,_agnxtout,0,__ZNKSt3__18time_putIcNS_19ostreambuf_iteratorIcNS_11char_traitsIcEEEEE6do_putES4_RNS_8ios_baseEcPK2tmcc,0,_tkgen_ellipse,0,__ZN10emscripten8internal14raw_destructorI8Agdesc_sEEvPT_,0,__ZNKSt3__120__time_get_c_storageIcE3__rEv
,0,_record_inside,0,_fig_begin_node,0,__ZNKSt3__114error_category10equivalentEiRKNS_15error_conditionE,0,_bothfunc,0,_nonefunc
,0,__ZNSt3__115basic_streambufIcNS_11char_traitsIcEEE6xsputnEPKci,0,_nop2_layout,0,_fig_end_edge,0,___cxx_global_array_dtor56,0,__ZNKSt3__15ctypeIwE10do_scan_isEtPKwS3_
,0,_freeMPair,0,__ZNKSt3__17codecvtIDic10_mbstate_tE6do_outERS1_PKDiS5_RS5_PcS7_RS7_,0,_aglstnode,0,_distX,0,__ZNKSt3__17codecvtIDic10_mbstate_tE5do_inERS1_PKcS5_RS5_PDiS7_RS7_
,0,_svg_ellipse,0,__ZTv0_n12_NSt3__113basic_ostreamIwNS_11char_traitsIwEEED0Ev,0,__ZNSt3__115basic_streambufIwNS_11char_traitsIwEEED1Ev,0,__ZN10emscripten8internal13getActualTypeI8Agnode_sEEPKNS0_7_TYPEIDEPT_,0,__ZN10__cxxabiv120__si_class_type_infoD0Ev
,0,__ZNKSt3__17collateIwE10do_compareEPKwS3_S3_S3_,0,_fig_polygon,0,___ZN10emscripten8internal7InvokerINSt3__112basic_stringIcNS2_11char_traitsIcEENS2_9allocatorIcEEEEJP8Agraph_sRKS8_EE6invokeEPFS8_SA_SC_ESA_PNS0_11BindingTypeIS8_E3$_0E_,0,_tkgen_bezier,0,_gvevent_motion
,0,__ZNSt3__115basic_streambufIcNS_11char_traitsIcEEE6xsgetnEPci,0,__ZNKSt3__17num_getIwNS_19istreambuf_iteratorIwNS_11char_traitsIwEEEEE6do_getES4_S4_RNS_8ios_baseERjRPv,0,_twopi_layout,0,_psgen_end_node,0,_dtlist
,0,__ZNKSt3__15ctypeIcE10do_tolowerEc,0,__ZNKSt3__110moneypunctIwLb1EE13do_neg_formatEv,0,_compFunction2,0,_tkgen_begin_job,0,__ZNKSt3__15ctypeIcE8do_widenEPKcS3_Pc
,0,__ZNSt3__17codecvtIwc10_mbstate_tED0Ev,0,_agerrorf,0,_ioputstr,0,__ZNKSt3__110moneypunctIwLb1EE16do_decimal_pointEv,0,_pic_polygon
,0,__ZNSt3__17codecvtIDsc10_mbstate_tED0Ev,0,_psgen_comment,0,__ZNKSt3__120__time_get_c_storageIcE7__weeksEv,0,_core_loadimage_ps,0,__ZNKSt3__17codecvtIcc10_mbstate_tE5do_inERS1_PKcS5_RS5_PcS7_RS7_
,0,__ZNKSt3__18numpunctIwE11do_truenameEv,0,_psgen_begin_edge,0,_epsf_inside,0,__ZTv0_n12_NSt3__113basic_istreamIcNS_11char_traitsIcEEED1Ev,0,_cmpItem3767
,0,__ZN10emscripten8internal14raw_destructorI8Agedge_sEEvPT_,0,_swap_ends_p4417,0,__ZNSt3__110__stdinbufIwE9underflowEv,0,__ZNSt3__18ios_base7failureD0Ev,0,_svg_textpara
,0,__ZNSt3__19money_getIcNS_19istreambuf_iteratorIcNS_11char_traitsIcEEEEED0Ev,0,__ZNSt3__18ios_base4InitD2Ev,0,__ZNKSt3__15ctypeIwE5do_isEtw,0,__ZNSt3__110moneypunctIwLb1EED0Ev,0,_pic_polyline
,0,__ZNKSt3__15ctypeIcE9do_narrowEPKcS3_cPc,0,__ZNKSt3__18time_getIcNS_19istreambuf_iteratorIcNS_11char_traitsIcEEEEE14do_get_weekdayES4_S4_RNS_8ios_baseERjP2tm,0,__ZNKSt3__17codecvtIDic10_mbstate_tE16do_always_noconvEv,0,__ZNKSt3__15ctypeIcE10do_toupperEPcPKc,0,___ZN10emscripten8internal7InvokerIiJP8Agnode_sRKNSt3__112basic_stringIcNS4_11char_traitsIcEENS4_9allocatorIcEEEESC_EE6invokeEPFiS3_SC_SC_ES3_PNS0_11BindingTypeISA_E3$_0ESJ__
,0,_psgen_end_job,0,_idclose,0,___cxx_global_array_dtor105,0,_memalloc,0,__ZNSt3__115basic_streambufIwNS_11char_traitsIwEEE7seekoffExNS_8ios_base7seekdirEj
,0,__ZNSt3__115basic_streambufIwNS_11char_traitsIwEEE6setbufEPwi,0,__ZNKSt3__18messagesIwE7do_openERKNS_12basic_stringIcNS_11char_traitsIcEENS_9allocatorIcEEEERKNS_6localeE,0,___ZN10emscripten8internal7InvokerIP8Agraph_sJS3_RKNSt3__112basic_stringIcNS4_11char_traitsIcEENS4_9allocatorIcEEEEEE6invokeEPFS3_S3_SC_ES3_PNS0_11BindingTypeISA_E3$_0E_,0,_comp_entities,0,__ZNSt3__17codecvtIDic10_mbstate_tED0Ev
,0,_vml_begin_graph,0,__ZN10emscripten8internal7InvokerIP8Agnode_sJP8Agraph_sPciEE6invokeEPFS3_S5_S6_iES5_S6_i,0,_ijcmpf,0,_subnode_search,0,_edgeidcmpf
,0,_pic_bezier,0,__ZNKSt3__110moneypunctIcLb1EE14do_curr_symbolEv,0,_fig_ellipse,0,__ZNSt3__16locale5__impD0Ev,0,__ZNK10__cxxabiv117__class_type_info27has_unambiguous_public_baseEPNS_19__dynamic_cast_infoEPvi
,0,__ZN10emscripten8internal7InvokerIP5GVC_sJEE6invokeEPFS3_vE,0,__ZNKSt3__119__iostream_category4nameEv,0,__ZNKSt3__110moneypunctIcLb0EE14do_frac_digitsEv,0,__ZNKSt3__18time_getIcNS_19istreambuf_iteratorIcNS_11char_traitsIcEEEEE16do_get_monthnameES4_S4_RNS_8ios_baseERjP2tm,0,__ZNKSt3__110moneypunctIwLb1EE11do_groupingEv
,0,__ZNSt3__19money_getIcNS_19istreambuf_iteratorIcNS_11char_traitsIcEEEEED1Ev,0,_xdot_polyline,0,_idalloc,0,_gvFreeLayout,0,_psgen_begin_page
,0,__ZNSt3__18time_putIcNS_19ostreambuf_iteratorIcNS_11char_traitsIcEEEEED1Ev,0,_agnxtnode,0,_svg_end_node,0,__ZNSt3__19money_getIwNS_19istreambuf_iteratorIwNS_11char_traitsIwEEEEED1Ev,0,__ZNSt8bad_castD0Ev
,0,__ZNKSt3__15ctypeIcE9do_narrowEcc,0,_agraphidcmpf,0,__ZNKSt3__17num_getIcNS_19istreambuf_iteratorIcNS_11char_traitsIcEEEEE6do_getES4_S4_RNS_8ios_baseERjRf,0,_fig_textpara,0,_mkItem3765
,0,__ZNSt3__112__do_nothingEPv,0,__ZNSt3__19money_putIwNS_19ostreambuf_iteratorIwNS_11char_traitsIwEEEEED1Ev,0,___cxx_global_array_dtor81,0,__ZNSt3__110moneypunctIcLb0EED0Ev,0,_newCell
,0,__ZNSt3__17num_putIwNS_19ostreambuf_iteratorIwNS_11char_traitsIwEEEEED0Ev,0,__ZNKSt3__17num_getIcNS_19istreambuf_iteratorIcNS_11char_traitsIcEEEEE6do_getES4_S4_RNS_8ios_baseERjRb,0,_vml_end_anchor,0,__ZNKSt3__17num_putIcNS_19ostreambuf_iteratorIcNS_11char_traitsIcEEEEE6do_putES4_RNS_8ios_baseEcPKv,0,_toggle_fit_cb
,0,__ZNKSt3__18numpunctIwE12do_falsenameEv,0,__ZNSt3__17collateIcED0Ev,0,__ZNKSt3__110moneypunctIwLb0EE13do_pos_formatEv,0,_psgen_begin_job,0,__ZNKSt3__110moneypunctIcLb1EE16do_negative_signEv
,0,_gvflush,0,__ZNSt3__111__stdoutbufIcED0Ev,0,_tkgen_begin_graph,0,_core_loadimage_pslib,0,__ZNSt3__16locale5facetD2Ev
,0,__ZTv0_n12_NSt3__113basic_istreamIwNS_11char_traitsIwEEED1Ev,0,__ZNSt3__112system_errorD0Ev,0,__ZNKSt3__19money_getIcNS_19istreambuf_iteratorIcNS_11char_traitsIcEEEEE6do_getES4_S4_bRNS_8ios_baseERjRe,0,__ZNKSt3__17codecvtIcc10_mbstate_tE9do_lengthERS1_PKcS5_j,0,__ZNSt3__110__stdinbufIwE5uflowEv
,0,__ZNKSt3__18numpunctIcE11do_truenameEv,0,_svg_bezier,0,__ZNKSt3__110moneypunctIcLb1EE13do_pos_formatEv,0,__ZNSt3__115basic_streambufIwNS_11char_traitsIwEEE7seekposENS_4fposI10_mbstate_tEEj,0,_poly_port
,0,_scomp,0,_pov_end_node,0,__Z9agnodesetP8Agnode_sRKNSt3__112basic_stringIcNS1_11char_traitsIcEENS1_9allocatorIcEEEES9_,0,_mkMPair,0,__ZNKSt3__19money_getIwNS_19istreambuf_iteratorIwNS_11char_traitsIwEEEEE6do_getES4_S4_bRNS_8ios_baseERjRe
,0,__ZNKSt3__17num_getIwNS_19istreambuf_iteratorIwNS_11char_traitsIwEEEEE6do_getES4_S4_RNS_8ios_baseERjS8_,0,_tkgen_textpara,0,__spline_edges,0,__ZNKSt3__18numpunctIcE16do_thousands_sepEv,0,__ZNSt3__110__stdinbufIwE5imbueERKNS_6localeE
,0,_memopen,0,_agidnode,0,__ZNSt3__19money_putIcNS_19ostreambuf_iteratorIcNS_11char_traitsIcEEEEED0Ev,0,__ZNSt3__113basic_istreamIcNS_11char_traitsIcEEED0Ev,0,_freePara
,0,_my_init_edge,0,__ZNSt3__115basic_streambufIcNS_11char_traitsIcEEE9showmanycEv,0,_intersectY,0,_intersectX,0,__ZNKSt3__18time_getIwNS_19istreambuf_iteratorIwNS_11char_traitsIwEEEEE11do_get_dateES4_S4_RNS_8ios_baseERjP2tm
,0,_free_clust,0,__ZNSt3__115basic_streambufIwNS_11char_traitsIwEEE8overflowEi,0,_pov_end_layer,0,_pov_comment,0,_psgen_begin_anchor
,0,__ZNSt3__18numpunctIwED0Ev,0,_idregister,0,__ZNKSt3__17num_getIwNS_19istreambuf_iteratorIwNS_11char_traitsIwEEEEE6do_getES4_S4_RNS_8ios_baseERjRx,0,_intersectX0,0,__ZNK10__cxxabiv119__pointer_type_info9can_catchEPKNS_16__shim_type_infoERPv
,0,__ZNSt3__115basic_streambufIcNS_11char_traitsIcEEE5imbueERKNS_6localeE,0,__ZNKSt3__15ctypeIwE10do_tolowerEw,0,__ZN10emscripten8internal7InvokerIiJP8Agraph_sEE6invokeEPFiS3_ES3_,0,__ZN10emscripten8internal7InvokerIiJP8Agraph_sP8Agnode_sEE6invokeEPFiS3_S5_ES3_S5_,0,_poly_free
,0,_fontcmpf,0,_gvevent_modify,0,_dot_end_graph,0,__ZNSt3__112basic_stringIwNS_11char_traitsIwEENS_9allocatorIwEEED1Ev,0,_agraphattr_init
,0,_tkgen_polygon,0,__Z17wrap_gvcBuildDateP5GVC_s,0,__ZNKSt3__17codecvtIwc10_mbstate_tE10do_unshiftERS1_PcS4_RS4_,0,_idcmpf,0,__ZNKSt3__17collateIcE10do_compareEPKcS3_S3_S3_
,0,__ZNKSt3__17codecvtIwc10_mbstate_tE16do_always_noconvEv,0,_idopen,0,_fig_resolve_color,0,__ZNKSt3__17codecvtIwc10_mbstate_tE6do_outERS1_PKwS5_RS5_PcS7_RS7_,0,_record_path
,0,__ZNSt3__110__stdinbufIwED1Ev,0,_Operator_matmul_apply,0,_psgen_polyline,0,__ZNKSt3__17collateIwE12do_transformEPKwS3_,0,__ZNKSt3__17num_putIcNS_19ostreambuf_iteratorIcNS_11char_traitsIcEEEEE6do_putES4_RNS_8ios_baseEcy
,0,__ZNKSt3__17num_putIcNS_19ostreambuf_iteratorIcNS_11char_traitsIcEEEEE6do_putES4_RNS_8ios_baseEcx,0,__ZNKSt3__17num_getIwNS_19istreambuf_iteratorIwNS_11char_traitsIwEEEEE6do_getES4_S4_RNS_8ios_baseERjRt,0,__ZNKSt3__17num_putIcNS_19ostreambuf_iteratorIcNS_11char_traitsIcEEEEE6do_putES4_RNS_8ios_baseEce,0,__ZNKSt3__17num_putIcNS_19ostreambuf_iteratorIcNS_11char_traitsIcEEEEE6do_putES4_RNS_8ios_baseEcd,0,_gvevent_read
,0,_record_gencode,0,_dttree,0,__ZNKSt3__17num_putIcNS_19ostreambuf_iteratorIcNS_11char_traitsIcEEEEE6do_putES4_RNS_8ios_baseEcl,0,__ZNSt8bad_castD2Ev,0,_vml_begin_anchor
,0,__ZN10__cxxabiv121__vmi_class_type_infoD0Ev,0,_scomp4544,0,_svg_end_edge,0,__ZNKSt3__110moneypunctIwLb1EE13do_pos_formatEv,0,_fdp_cleanup
,0,__ZNKSt3__110moneypunctIcLb1EE14do_frac_digitsEv,0,_pov_polyline,0,__ZNKSt3__17codecvtIDic10_mbstate_tE10do_unshiftERS1_PcS4_RS4_,0,__ZNKSt3__19money_getIcNS_19istreambuf_iteratorIcNS_11char_traitsIcEEEEE6do_getES4_S4_bRNS_8ios_baseERjRNS_12basic_stringIcS3_NS_9allocatorIcEEEE,0,__ZNKSt3__15ctypeIwE10do_toupperEPwPKw
,0,__ZTv0_n12_NSt3__113basic_ostreamIcNS_11char_traitsIcEEED0Ev,0,_agnxtedge,0,__ZNSt3__110__stdinbufIcE9underflowEv,0,__ZNKSt3__17num_getIwNS_19istreambuf_iteratorIwNS_11char_traitsIwEEEEE6do_getES4_S4_RNS_8ios_baseERjRl,0,__ZNKSt3__114error_category23default_error_conditionEi
,0,_free_item,0,_vml_comment,0,_up_cb,0,__ZNKSt3__17codecvtIcc10_mbstate_tE13do_max_lengthEv,0,__ZNK10__cxxabiv117__class_type_info16search_above_dstEPNS_19__dynamic_cast_infoEPKvS4_ib
,0,__ZN10emscripten8internal7InvokerIiJP8Agraph_sP8Agedge_sEE6invokeEPFiS3_S5_ES3_S5_,0,_svg_begin_anchor,0,_point_gencode,0,__ZNKSt3__17codecvtIcc10_mbstate_tE16do_always_noconvEv,0,__ZNKSt3__18messagesIwE8do_closeEi
,0,_neato_cleanup,0,_dot_layout,0,_arrow_type_crow,0,__ZNSt3__110__stdinbufIcE5imbueERKNS_6localeE,0,__ZNSt3__112system_errorD2Ev
,0,__ZNKSt9bad_alloc4whatEv,0,__ZNSt3__15ctypeIwED0Ev,0,__ZNKSt3__110moneypunctIwLb0EE11do_groupingEv,0,__ZNK10__cxxabiv123__fundamental_type_info9can_catchEPKNS_16__shim_type_infoERPv,0,__ZNSt3__115basic_streambufIwNS_11char_traitsIwEEE9showmanycEv
,0,__ZNKSt3__110moneypunctIcLb0EE16do_decimal_pointEv,0,__Z11wrap_agedgeP8Agraph_sP8Agnode_sS2_RKNSt3__112basic_stringIcNS3_11char_traitsIcEENS3_9allocatorIcEEEEi,0,_pov_bezier,0,__ZNKSt3__17num_getIwNS_19istreambuf_iteratorIwNS_11char_traitsIwEEEEE6do_getES4_S4_RNS_8ios_baseERjRy,0,__Z15wrap_gvcVersionP5GVC_s
,0,_free_subnode,0,_gvevent_refresh,0,_vml_bezier,0,__ZNKSt3__120__time_get_c_storageIcE8__monthsEv,0,_idprint
,0,__ZNKSt3__17num_getIcNS_19istreambuf_iteratorIcNS_11char_traitsIcEEEEE6do_getES4_S4_RNS_8ios_baseERjRPv,0,_left_cb,0,_core_loadimage_fig,0,__ZNKSt3__17num_getIwNS_19istreambuf_iteratorIwNS_11char_traitsIwEEEEE6do_getES4_S4_RNS_8ios_baseERjRm,0,_agclose
,0,_tkgen_comment,0,_free,0,___ZN10emscripten8internal7InvokerIP8Agedge_sJP8Agraph_sP8Agnode_sS7_RKNSt3__112basic_stringIcNS8_11char_traitsIcEENS8_9allocatorIcEEEEiEE6invokeEPFS3_S5_S7_S7_SG_iES5_S7_S7_PNS0_11BindingTypeISE_E3$_0Ei_,0,__ZNKSt3__17num_getIwNS_19istreambuf_iteratorIwNS_11char_traitsIwEEEEE6do_getES4_S4_RNS_8ios_baseERjRb,0,__ZNKSt3__17num_getIwNS_19istreambuf_iteratorIwNS_11char_traitsIwEEEEE6do_getES4_S4_RNS_8ios_baseERjRe
,0,__ZNKSt3__17num_getIwNS_19istreambuf_iteratorIwNS_11char_traitsIwEEEEE6do_getES4_S4_RNS_8ios_baseERjRd,0,_epsf_free,0,__ZNKSt3__17num_getIwNS_19istreambuf_iteratorIwNS_11char_traitsIwEEEEE6do_getES4_S4_RNS_8ios_baseERjRf,0,__ZNKSt3__17codecvtIcc10_mbstate_tE11do_encodingEv,0,_spline_merge
,0,_pic_ellipse,0,_comp_real,0,_circo_cleanup,0,__ZNKSt3__114error_category10equivalentERKNS_10error_codeEi,0,__ZNKSt3__110moneypunctIcLb0EE13do_neg_formatEv
,0,__ZNKSt11logic_error4whatEv,0,__ZNKSt3__119__iostream_category7messageEi,0,__ZNKSt3__17codecvtIDsc10_mbstate_tE13do_max_lengthEv,0,__ZNSt11logic_errorD0Ev,0,__ZNKSt3__110moneypunctIcLb0EE16do_thousands_sepEv
,0,_agnxtin,0,__ZNKSt3__110moneypunctIwLb1EE14do_curr_symbolEv,0,__ZNKSt3__110moneypunctIcLb0EE13do_pos_formatEv,0,_core_loadimage_xdot,0,_pov_begin_page
,0,__ZNSt3__113basic_ostreamIwNS_11char_traitsIwEEED0Ev,0,__ZNKSt3__110moneypunctIwLb0EE16do_decimal_pointEv,0,__ZNKSt3__17codecvtIDsc10_mbstate_tE16do_always_noconvEv,0,__ZNKSt3__17collateIcE12do_transformEPKcS3_,0,__ZNSt3__115basic_streambufIwNS_11char_traitsIwEEE6xsgetnEPwi
,0,_compFunction,0,__Z14gvRenderStringP5GVC_sP8Agraph_sRKNSt3__112basic_stringIcNS3_11char_traitsIcEENS3_9allocatorIcEEEE,0,__ZNKSt3__110moneypunctIwLb0EE14do_frac_digitsEv,0,__ZN10emscripten8internal7InvokerIiJP8Agraph_sP8Agnode_siiEE6invokeEPFiS3_S5_iiES3_S5_ii,0,__ZN10emscripten8internal14raw_destructorI8Agnode_sEEvPT_
,0,__ZNKSt3__110moneypunctIwLb0EE16do_thousands_sepEv,0,_svg_begin_graph,0,__ZNKSt3__15ctypeIcE10do_tolowerEPcPKc,0,_arrow_type_tee,0,__ZNKSt3__120__time_get_c_storageIcE7__am_pmEv
,0,__ZNKSt3__110moneypunctIcLb0EE14do_curr_symbolEv,0,__ZNKSt3__15ctypeIwE8do_widenEPKcS3_Pw,0,__ZNKSt3__110moneypunctIwLb1EE16do_thousands_sepEv,0,_gvevent_delete,0,__ZNK10__cxxabiv121__vmi_class_type_info16search_below_dstEPNS_19__dynamic_cast_infoEPKvib
,0,___cxx_global_array_dtor53,0,___ZN10emscripten8internal7InvokerIP8Agraph_sJRKNSt3__112basic_stringIcNS4_11char_traitsIcEENS4_9allocatorIcEEEE8Agdesc_sEE6invokeEPFS3_SC_SD_EPNS0_11BindingTypeISA_E3$_0EPSD__,0,_namecmpf,0,_newItem,0,_ps_freeimage
,0,__ZNSt3__18ios_baseD2Ev,0,__ZNSt3__113basic_ostreamIcNS_11char_traitsIcEEED1Ev,0,_vml_polygon,0,__ZNSt3__110__stdinbufIcED0Ev,0,_agfstout
,0,__ZNKSt3__17num_putIcNS_19ostreambuf_iteratorIcNS_11char_traitsIcEEEEE6do_putES4_RNS_8ios_baseEcb,0,__ZN10emscripten8internal7InvokerIP8Agedge_sJP8Agraph_sS3_P8Agnode_sEE6invokeEPFS3_S5_S3_S7_ES5_S3_S7_,0,_free_string_entry,0,__ZNKSt3__18time_getIwNS_19istreambuf_iteratorIwNS_11char_traitsIwEEEEE11do_get_yearES4_S4_RNS_8ios_baseERjP2tm,0,_freesym
,0,_map_end_page,0,__ZNSt3__110moneypunctIwLb1EED1Ev,0,_idmap,0,__ZNKSt3__110moneypunctIwLb0EE14do_curr_symbolEv,0,_agfstedge
,0,__ZNKSt3__18messagesIcE7do_openERKNS_12basic_stringIcNS_11char_traitsIcEENS_9allocatorIcEEEERKNS_6localeE,0,_addattr,0,_cmpitem,0,__Z9agedgesetP8Agedge_sRKNSt3__112basic_stringIcNS1_11char_traitsIcEENS1_9allocatorIcEEEES9_,0,__ZNSt3__18time_getIcNS_19istreambuf_iteratorIcNS_11char_traitsIcEEEEED0Ev
,0,__ZNSt3__110moneypunctIcLb1EED1Ev,0,_psgen_begin_cluster,0,__ZNSt3__111__stdoutbufIwED0Ev,0,__ZNSt3__115basic_streambufIcNS_11char_traitsIcEEE7seekoffExNS_8ios_base7seekdirEj,0,__ZNKSt3__120__time_get_c_storageIcE3__cEv
,0,__ZNSt3__17codecvtIwc10_mbstate_tED2Ev,0,__ZNSt3__115basic_streambufIcNS_11char_traitsIcEEE6setbufEPci,0,_svg_end_graph,0,__ZNKSt3__110moneypunctIwLb0EE16do_positive_signEv,0,_distY
,0,__ZN10emscripten8internal7InvokerIP8Agnode_sJP8Agraph_sEE6invokeEPFS3_S5_ES5_,0,__ZNKSt3__17num_getIcNS_19istreambuf_iteratorIcNS_11char_traitsIcEEEEE6do_getES4_S4_RNS_8ios_baseERjS8_,0,_nop1_layout,0,__ZNKSt3__17codecvtIDic10_mbstate_tE13do_max_lengthEv,0,__ZNKSt3__120__time_get_c_storageIwE3__XEv
,0,_forfunc,0,_ioflush,0,__Z13wrap_agconcatP8Agraph_sRKNSt3__112basic_stringIcNS1_11char_traitsIcEENS1_9allocatorIcEEEE,0,__ZNKSt3__18time_getIcNS_19istreambuf_iteratorIcNS_11char_traitsIcEEEEE11do_get_dateES4_S4_RNS_8ios_baseERjP2tm,0,__ZNKSt3__110moneypunctIcLb1EE16do_positive_signEv
,0,_gvevent_render,0,__ZTv0_n12_NSt3__113basic_istreamIwNS_11char_traitsIwEEED0Ev,0,__ZNSt3__115basic_streambufIwNS_11char_traitsIwEEED0Ev,0,_usershape_close,0,_freeItem
,0,__ZNKSt3__18numpunctIwE11do_groupingEv,0,__ZNSt3__115basic_streambufIcNS_11char_traitsIcEEE9pbackfailEi,0,__ZNSt3__113basic_ostreamIcNS_11char_traitsIcEEED0Ev,0,_quit_cb,0,_makePoly
,0,__ZN10emscripten8internal7InvokerIP8Agnode_sJP8Agraph_sS3_EE6invokeEPFS3_S5_S3_ES5_S3_,0,__ZNSt3__111__stdoutbufIwE8overflowEi,0,__ZNKSt3__17num_getIcNS_19istreambuf_iteratorIcNS_11char_traitsIcEEEEE6do_getES4_S4_RNS_8ios_baseERjRy,0,__ZNKSt3__17num_getIcNS_19istreambuf_iteratorIcNS_11char_traitsIcEEEEE6do_getES4_S4_RNS_8ios_baseERjRx,0,__ZNKSt3__17num_getIcNS_19istreambuf_iteratorIcNS_11char_traitsIcEEEEE6do_getES4_S4_RNS_8ios_baseERjRt
,0,__ZNKSt3__17num_getIcNS_19istreambuf_iteratorIcNS_11char_traitsIcEEEEE6do_getES4_S4_RNS_8ios_baseERjRm,0,__ZNKSt3__17num_getIcNS_19istreambuf_iteratorIcNS_11char_traitsIcEEEEE6do_getES4_S4_RNS_8ios_baseERjRl,0,__ZNKSt3__19money_putIcNS_19ostreambuf_iteratorIcNS_11char_traitsIcEEEEE6do_putES4_bRNS_8ios_baseEce,0,_swap_ends_p,0,__ZNKSt3__17num_getIcNS_19istreambuf_iteratorIcNS_11char_traitsIcEEEEE6do_getES4_S4_RNS_8ios_baseERjRe
,0,__ZNKSt3__17num_getIcNS_19istreambuf_iteratorIcNS_11char_traitsIcEEEEE6do_getES4_S4_RNS_8ios_baseERjRd,0,__ZNSt3__116__narrow_to_utf8ILj32EED0Ev,0,_gvFreeContext,0,_poly_gencode,0,__ZN10emscripten8internal13getActualTypeI5GVC_sEEPKNS0_7_TYPEIDEPT_
,0,___cxx_global_array_dtor,0,_pov_ellipse,0,_poly_inside,0,__ZNKSt3__19money_putIwNS_19ostreambuf_iteratorIwNS_11char_traitsIwEEEEE6do_putES4_bRNS_8ios_baseEwe,0,_agdelnode
,0,_gvevent_button_release,0,__ZNSt3__18time_getIcNS_19istreambuf_iteratorIcNS_11char_traitsIcEEEEED1Ev,0,_agdictobjmem,0,__ZN10__cxxabiv117__class_type_infoD0Ev,0,__ZNSt3__18messagesIwED1Ev
,0,__ZNSt3__111__stdoutbufIwED1Ev,0,__Z11wrap_agopenRKNSt3__112basic_stringIcNS_11char_traitsIcEENS_9allocatorIcEEEE8Agdesc_s,0,__ZNKSt3__19money_getIwNS_19istreambuf_iteratorIwNS_11char_traitsIwEEEEE6do_getES4_S4_bRNS_8ios_baseERjRNS_12basic_stringIwS3_NS_9allocatorIwEEEE,0,_pov_end_graph,0,_edgecmp
,0,__ZN10__cxxabiv116__shim_type_infoD2Ev,0,__ZNKSt3__15ctypeIwE11do_scan_notEtPKwS3_,0,__ZNSt3__18time_putIwNS_19ostreambuf_iteratorIwNS_11char_traitsIwEEEEED1Ev,0,_inside,0,_svg_polygon
,0,_svg_end_anchor,0,__ZNKSt3__17codecvtIDsc10_mbstate_tE6do_outERS1_PKDsS5_RS5_PcS7_RS7_,0,_intersectY0,0,__ZNKSt3__18messagesIwE6do_getEiiiRKNS_12basic_stringIwNS_11char_traitsIwEENS_9allocatorIwEEEE,0,__ZNKSt3__110moneypunctIcLb1EE11do_groupingEv
,0,___getTypeName,0,_svg_comment,0,__ZNKSt3__17codecvtIDsc10_mbstate_tE11do_encodingEv,0,__ZNSt3__115basic_streambufIwNS_11char_traitsIwEEE9pbackfailEi,0,_freeItem3756
,0,_psgen_begin_layer,0,_idfree,0,__ZNSt3__115basic_streambufIwNS_11char_traitsIwEEE6xsputnEPKwi,0,_revfunc,0,__ZNSt3__112basic_stringIcNS_11char_traitsIcEENS_9allocatorIcEEED1Ev
,0,__ZNSt3__15ctypeIcED2Ev,0,_patchwork_layout,0,_Operator_uniform_stress_matmul_apply,0,__ZNSt13runtime_errorD0Ev,0,__ZNSt3__113basic_istreamIwNS_11char_traitsIwEEED0Ev,0,__ZN10emscripten8internal7InvokerIP8Agedge_sJP8Agraph_sS3_EE6invokeEPFS3_S5_S3_ES5_S3_,0,___cxx_global_array_dtor120];
// EMSCRIPTEN_START_FUNCS
function _dispose_chunk(r1,r2){var r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15,r16,r17,r18,r19,r20,r21,r22,r23,r24,r25,r26,r27,r28,r29,r30,r31,r32,r33,r34,r35,r36,r37,r38,r39,r40,r41,r42;r3=r2>>2;r4=0;r5=r1,r6=r5>>2;r7=r5+r2|0;r8=r7;r9=HEAP32[r1+4>>2];L59006:do{if((r9&1|0)==0){r10=HEAP32[r1>>2];if((r9&3|0)==0){return}r11=r5+ -r10|0;r12=r11;r13=r10+r2|0;r14=HEAP32[43906];if(r11>>>0<r14>>>0){_abort()}if((r12|0)==(HEAP32[43907]|0)){r15=(r2+(r5+4)|0)>>2;if((HEAP32[r15]&3|0)!=3){r16=r12,r17=r16>>2;r18=r13;break}HEAP32[43904]=r13;HEAP32[r15]=HEAP32[r15]&-2;HEAP32[(4-r10>>2)+r6]=r13|1;HEAP32[r7>>2]=r13;return}r15=r10>>>3;if(r10>>>0<256){r19=HEAP32[(8-r10>>2)+r6];r20=HEAP32[(12-r10>>2)+r6];r21=(r15<<3)+175648|0;do{if((r19|0)!=(r21|0)){if(r19>>>0<r14>>>0){_abort()}if((HEAP32[r19+12>>2]|0)==(r12|0)){break}_abort()}}while(0);if((r20|0)==(r19|0)){HEAP32[43902]=HEAP32[43902]&(1<<r15^-1);r16=r12,r17=r16>>2;r18=r13;break}do{if((r20|0)==(r21|0)){r22=r20+8|0}else{if(r20>>>0<r14>>>0){_abort()}r23=r20+8|0;if((HEAP32[r23>>2]|0)==(r12|0)){r22=r23;break}_abort()}}while(0);HEAP32[r19+12>>2]=r20;HEAP32[r22>>2]=r19;r16=r12,r17=r16>>2;r18=r13;break}r21=r11;r15=HEAP32[(24-r10>>2)+r6];r23=HEAP32[(12-r10>>2)+r6];do{if((r23|0)==(r21|0)){r24=16-r10|0;r25=r24+(r5+4)|0;r26=HEAP32[r25>>2];if((r26|0)==0){r27=r5+r24|0;r24=HEAP32[r27>>2];if((r24|0)==0){r28=0,r29=r28>>2;break}else{r30=r24;r31=r27}}else{r30=r26;r31=r25}while(1){r25=r30+20|0;r26=HEAP32[r25>>2];if((r26|0)!=0){r30=r26;r31=r25;continue}r25=r30+16|0;r26=HEAP32[r25>>2];if((r26|0)==0){break}else{r30=r26;r31=r25}}if(r31>>>0<r14>>>0){_abort()}else{HEAP32[r31>>2]=0;r28=r30,r29=r28>>2;break}}else{r25=HEAP32[(8-r10>>2)+r6];if(r25>>>0<r14>>>0){_abort()}r26=r25+12|0;if((HEAP32[r26>>2]|0)!=(r21|0)){_abort()}r27=r23+8|0;if((HEAP32[r27>>2]|0)==(r21|0)){HEAP32[r26>>2]=r23;HEAP32[r27>>2]=r25;r28=r23,r29=r28>>2;break}else{_abort()}}}while(0);if((r15|0)==0){r16=r12,r17=r16>>2;r18=r13;break}r23=r5+(28-r10)|0;r14=(HEAP32[r23>>2]<<2)+175912|0;do{if((r21|0)==(HEAP32[r14>>2]|0)){HEAP32[r14>>2]=r28;if((r28|0)!=0){break}HEAP32[43903]=HEAP32[43903]&(1<<HEAP32[r23>>2]^-1);r16=r12,r17=r16>>2;r18=r13;break L59006}else{if(r15>>>0<HEAP32[43906]>>>0){_abort()}r11=r15+16|0;if((HEAP32[r11>>2]|0)==(r21|0)){HEAP32[r11>>2]=r28}else{HEAP32[r15+20>>2]=r28}if((r28|0)==0){r16=r12,r17=r16>>2;r18=r13;break L59006}}}while(0);if(r28>>>0<HEAP32[43906]>>>0){_abort()}HEAP32[r29+6]=r15;r21=16-r10|0;r23=HEAP32[(r21>>2)+r6];do{if((r23|0)!=0){if(r23>>>0<HEAP32[43906]>>>0){_abort()}else{HEAP32[r29+4]=r23;HEAP32[r23+24>>2]=r28;break}}}while(0);r23=HEAP32[(r21+4>>2)+r6];if((r23|0)==0){r16=r12,r17=r16>>2;r18=r13;break}if(r23>>>0<HEAP32[43906]>>>0){_abort()}else{HEAP32[r29+5]=r23;HEAP32[r23+24>>2]=r28;r16=r12,r17=r16>>2;r18=r13;break}}else{r16=r1,r17=r16>>2;r18=r2}}while(0);r1=HEAP32[43906];if(r7>>>0<r1>>>0){_abort()}r28=r2+(r5+4)|0;r29=HEAP32[r28>>2];do{if((r29&2|0)==0){if((r8|0)==(HEAP32[43908]|0)){r30=HEAP32[43905]+r18|0;HEAP32[43905]=r30;HEAP32[43908]=r16;HEAP32[r17+1]=r30|1;if((r16|0)!=(HEAP32[43907]|0)){return}HEAP32[43907]=0;HEAP32[43904]=0;return}if((r8|0)==(HEAP32[43907]|0)){r30=HEAP32[43904]+r18|0;HEAP32[43904]=r30;HEAP32[43907]=r16;HEAP32[r17+1]=r30|1;HEAP32[(r30>>2)+r17]=r30;return}r30=(r29&-8)+r18|0;r31=r29>>>3;L59105:do{if(r29>>>0<256){r22=HEAP32[r3+(r6+2)];r9=HEAP32[r3+(r6+3)];r23=(r31<<3)+175648|0;do{if((r22|0)!=(r23|0)){if(r22>>>0<r1>>>0){_abort()}if((HEAP32[r22+12>>2]|0)==(r8|0)){break}_abort()}}while(0);if((r9|0)==(r22|0)){HEAP32[43902]=HEAP32[43902]&(1<<r31^-1);break}do{if((r9|0)==(r23|0)){r32=r9+8|0}else{if(r9>>>0<r1>>>0){_abort()}r10=r9+8|0;if((HEAP32[r10>>2]|0)==(r8|0)){r32=r10;break}_abort()}}while(0);HEAP32[r22+12>>2]=r9;HEAP32[r32>>2]=r22}else{r23=r7;r10=HEAP32[r3+(r6+6)];r15=HEAP32[r3+(r6+3)];do{if((r15|0)==(r23|0)){r14=r2+(r5+20)|0;r11=HEAP32[r14>>2];if((r11|0)==0){r19=r2+(r5+16)|0;r20=HEAP32[r19>>2];if((r20|0)==0){r33=0,r34=r33>>2;break}else{r35=r20;r36=r19}}else{r35=r11;r36=r14}while(1){r14=r35+20|0;r11=HEAP32[r14>>2];if((r11|0)!=0){r35=r11;r36=r14;continue}r14=r35+16|0;r11=HEAP32[r14>>2];if((r11|0)==0){break}else{r35=r11;r36=r14}}if(r36>>>0<r1>>>0){_abort()}else{HEAP32[r36>>2]=0;r33=r35,r34=r33>>2;break}}else{r14=HEAP32[r3+(r6+2)];if(r14>>>0<r1>>>0){_abort()}r11=r14+12|0;if((HEAP32[r11>>2]|0)!=(r23|0)){_abort()}r19=r15+8|0;if((HEAP32[r19>>2]|0)==(r23|0)){HEAP32[r11>>2]=r15;HEAP32[r19>>2]=r14;r33=r15,r34=r33>>2;break}else{_abort()}}}while(0);if((r10|0)==0){break}r15=r2+(r5+28)|0;r22=(HEAP32[r15>>2]<<2)+175912|0;do{if((r23|0)==(HEAP32[r22>>2]|0)){HEAP32[r22>>2]=r33;if((r33|0)!=0){break}HEAP32[43903]=HEAP32[43903]&(1<<HEAP32[r15>>2]^-1);break L59105}else{if(r10>>>0<HEAP32[43906]>>>0){_abort()}r9=r10+16|0;if((HEAP32[r9>>2]|0)==(r23|0)){HEAP32[r9>>2]=r33}else{HEAP32[r10+20>>2]=r33}if((r33|0)==0){break L59105}}}while(0);if(r33>>>0<HEAP32[43906]>>>0){_abort()}HEAP32[r34+6]=r10;r23=HEAP32[r3+(r6+4)];do{if((r23|0)!=0){if(r23>>>0<HEAP32[43906]>>>0){_abort()}else{HEAP32[r34+4]=r23;HEAP32[r23+24>>2]=r33;break}}}while(0);r23=HEAP32[r3+(r6+5)];if((r23|0)==0){break}if(r23>>>0<HEAP32[43906]>>>0){_abort()}else{HEAP32[r34+5]=r23;HEAP32[r23+24>>2]=r33;break}}}while(0);HEAP32[r17+1]=r30|1;HEAP32[(r30>>2)+r17]=r30;if((r16|0)!=(HEAP32[43907]|0)){r37=r30;break}HEAP32[43904]=r30;return}else{HEAP32[r28>>2]=r29&-2;HEAP32[r17+1]=r18|1;HEAP32[(r18>>2)+r17]=r18;r37=r18}}while(0);r18=r37>>>3;if(r37>>>0<256){r29=r18<<1;r28=(r29<<2)+175648|0;r33=HEAP32[43902];r34=1<<r18;do{if((r33&r34|0)==0){HEAP32[43902]=r33|r34;r38=r28;r39=(r29+2<<2)+175648|0}else{r18=(r29+2<<2)+175648|0;r6=HEAP32[r18>>2];if(r6>>>0>=HEAP32[43906]>>>0){r38=r6;r39=r18;break}_abort()}}while(0);HEAP32[r39>>2]=r16;HEAP32[r38+12>>2]=r16;HEAP32[r17+2]=r38;HEAP32[r17+3]=r28;return}r28=r16;r38=r37>>>8;do{if((r38|0)==0){r40=0}else{if(r37>>>0>16777215){r40=31;break}r39=(r38+1048320|0)>>>16&8;r29=r38<<r39;r34=(r29+520192|0)>>>16&4;r33=r29<<r34;r29=(r33+245760|0)>>>16&2;r18=14-(r34|r39|r29)+(r33<<r29>>>15)|0;r40=r37>>>((r18+7|0)>>>0)&1|r18<<1}}while(0);r38=(r40<<2)+175912|0;HEAP32[r17+7]=r40;HEAP32[r17+5]=0;HEAP32[r17+4]=0;r18=HEAP32[43903];r29=1<<r40;if((r18&r29|0)==0){HEAP32[43903]=r18|r29;HEAP32[r38>>2]=r28;HEAP32[r17+6]=r38;HEAP32[r17+3]=r16;HEAP32[r17+2]=r16;return}if((r40|0)==31){r41=0}else{r41=25-(r40>>>1)|0}r40=r37<<r41;r41=HEAP32[r38>>2];while(1){if((HEAP32[r41+4>>2]&-8|0)==(r37|0)){break}r42=(r40>>>31<<2)+r41+16|0;r38=HEAP32[r42>>2];if((r38|0)==0){r4=44625;break}else{r40=r40<<1;r41=r38}}if(r4==44625){if(r42>>>0<HEAP32[43906]>>>0){_abort()}HEAP32[r42>>2]=r28;HEAP32[r17+6]=r41;HEAP32[r17+3]=r16;HEAP32[r17+2]=r16;return}r16=r41+8|0;r42=HEAP32[r16>>2];r4=HEAP32[43906];if(r41>>>0<r4>>>0){_abort()}if(r42>>>0<r4>>>0){_abort()}HEAP32[r42+12>>2]=r28;HEAP32[r16>>2]=r28;HEAP32[r17+2]=r42;HEAP32[r17+3]=r41;HEAP32[r17+6]=0;return}function __ZSt17__throw_bad_allocv(){var r1;r1=___cxa_allocate_exception(4);HEAP32[r1>>2]=180544;___cxa_throw(r1,187048,150)}function _strtod(r1,r2){var r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15,r16,r17,r18,r19,r20,r21,r22,r23,r24,r25,r26,r27,r28,r29,r30,r31,r32,r33,r34,r35,r36,r37,r38,r39,r40,r41;r3=0;r4=r1;while(1){r5=r4+1|0;if((_isspace(HEAP8[r4]|0)|0)==0){break}else{r4=r5}}r6=HEAP8[r4];if(r6<<24>>24==45){r7=r5;r8=1}else if(r6<<24>>24==43){r7=r5;r8=0}else{r7=r4;r8=0}r4=-1;r5=0;r6=r7;while(1){r9=HEAP8[r6];if(((r9<<24>>24)-48|0)>>>0<10){r10=r4}else{if(r9<<24>>24!=46|(r4|0)>-1){break}else{r10=r5}}r4=r10;r5=r5+1|0;r6=r6+1|0}r10=r6+ -r5|0;r7=(r4|0)<0;r11=((r7^1)<<31>>31)+r5|0;r12=(r11|0)>18;r13=(r12?-18:-r11|0)+(r7?r5:r4)|0;r4=r12?18:r11;do{if((r4|0)==0){r14=r1;r15=0}else{if((r4|0)>9){r11=r10;r12=r4;r5=0;while(1){r7=HEAP8[r11];r16=r11+1|0;if(r7<<24>>24==46){r17=HEAP8[r16];r18=r11+2|0}else{r17=r7;r18=r16}r19=(r17<<24>>24)+((r5*10&-1)-48)|0;r16=r12-1|0;if((r16|0)>9){r11=r18;r12=r16;r5=r19}else{break}}r20=(r19|0)*1e9;r21=9;r22=r18;r3=44674}else{if((r4|0)>0){r20=0;r21=r4;r22=r10;r3=44674}else{r23=0;r24=0}}if(r3==44674){r5=r22;r12=r21;r11=0;while(1){r16=HEAP8[r5];r7=r5+1|0;if(r16<<24>>24==46){r25=HEAP8[r7];r26=r5+2|0}else{r25=r16;r26=r7}r27=(r25<<24>>24)+((r11*10&-1)-48)|0;r7=r12-1|0;if((r7|0)>0){r5=r26;r12=r7;r11=r27}else{break}}r23=r27|0;r24=r20}r11=r24+r23;do{if(r9<<24>>24==69|r9<<24>>24==101){r12=r6+1|0;r5=HEAP8[r12];if(r5<<24>>24==45){r28=r6+2|0;r29=1}else if(r5<<24>>24==43){r28=r6+2|0;r29=0}else{r28=r12;r29=0}r12=HEAP8[r28];if(((r12<<24>>24)-48|0)>>>0<10){r30=r28;r31=0;r32=r12}else{r33=0;r34=r28;r35=r29;break}while(1){r12=(r32<<24>>24)+((r31*10&-1)-48)|0;r5=r30+1|0;r7=HEAP8[r5];if(((r7<<24>>24)-48|0)>>>0<10){r30=r5;r31=r12;r32=r7}else{r33=r12;r34=r5;r35=r29;break}}}else{r33=0;r34=r6;r35=0}}while(0);r5=r13+((r35|0)==0?r33:-r33|0)|0;r12=(r5|0)<0?-r5|0:r5;if((r12|0)>511){r7=___errno_location();HEAP32[r7>>2]=34;r36=1;r37=5400;r38=511;r3=44691}else{if((r12|0)==0){r39=1}else{r36=1;r37=5400;r38=r12;r3=44691}}if(r3==44691){while(1){r3=0;if((r38&1|0)==0){r40=r36}else{r40=r36*HEAPF64[r37>>3]}r12=r38>>1;if((r12|0)==0){r39=r40;break}else{r36=r40;r37=r37+8|0;r38=r12;r3=44691}}}if((r5|0)>-1){r14=r34;r15=r11*r39;break}else{r14=r34;r15=r11/r39;break}}}while(0);if((r2|0)!=0){HEAP32[r2>>2]=r14}if((r8|0)==0){r41=r15;return r41}r41=-r15;return r41}
// EMSCRIPTEN_END_FUNCS
Module["_main"] = _main;
Module["_malloc"] = _malloc;
Module["_realloc"] = _realloc;
// Warning: printing of i64 values may be slightly rounded! No deep i64 math used, so precise i64 code not included
var i64Math = null;
// === Auto-generated postamble setup entry stuff ===
Module['callMain'] = function callMain(args) {
  assert(runDependencies == 0, 'cannot call main when async dependencies remain! (listen on __ATMAIN__)');
  assert(!Module['preRun'] || Module['preRun'].length == 0, 'cannot call main when preRun functions remain to be called');
  args = args || [];
  ensureInitRuntime();
  var argc = args.length+1;
  function pad() {
    for (var i = 0; i < 4-1; i++) {
      argv.push(0);
    }
  }
  var argv = [allocate(intArrayFromString("/bin/this.program"), 'i8', ALLOC_NORMAL) ];
  pad();
  for (var i = 0; i < argc-1; i = i + 1) {
    argv.push(allocate(intArrayFromString(args[i]), 'i8', ALLOC_NORMAL));
    pad();
  }
  argv.push(0);
  argv = allocate(argv, 'i32', ALLOC_NORMAL);
  var ret;
  var initialStackTop = STACKTOP;
  try {
    ret = Module['_main'](argc, argv, 0);
  }
  catch(e) {
    if (e.name == 'ExitStatus') {
      return e.status;
    } else if (e == 'SimulateInfiniteLoop') {
      Module['noExitRuntime'] = true;
    } else {
      throw e;
    }
  } finally {
    STACKTOP = initialStackTop;
  }
  return ret;
}
function run(args) {
  args = args || Module['arguments'];
  if (runDependencies > 0) {
    Module.printErr('run() called, but dependencies remain, so not running');
    return 0;
  }
  if (Module['preRun']) {
    if (typeof Module['preRun'] == 'function') Module['preRun'] = [Module['preRun']];
    var toRun = Module['preRun'];
    Module['preRun'] = [];
    for (var i = toRun.length-1; i >= 0; i--) {
      toRun[i]();
    }
    if (runDependencies > 0) {
      // a preRun added a dependency, run will be called later
      return 0;
    }
  }
  function doRun() {
    ensureInitRuntime();
    preMain();
    var ret = 0;
    calledRun = true;
    if (Module['_main'] && shouldRunNow) {
      ret = Module['callMain'](args);
      if (!Module['noExitRuntime']) {
        exitRuntime();
      }
    }
    if (Module['postRun']) {
      if (typeof Module['postRun'] == 'function') Module['postRun'] = [Module['postRun']];
      while (Module['postRun'].length > 0) {
        Module['postRun'].pop()();
      }
    }
    return ret;
  }
  if (Module['setStatus']) {
    Module['setStatus']('Running...');
    setTimeout(function() {
      setTimeout(function() {
        Module['setStatus']('');
      }, 1);
      if (!ABORT) doRun();
    }, 1);
    return 0;
  } else {
    return doRun();
  }
}
Module['run'] = Module.run = run;
// {{PRE_RUN_ADDITIONS}}
/*global Module*/
/*global _malloc, _free, _memcpy*/
/*global FUNCTION_TABLE, HEAP8, HEAPU8, HEAP16, HEAPU16, HEAP32, HEAPU32*/
/*global readLatin1String*/
/*global __emval_register, _emval_handle_array, __emval_decref*/
/*global ___getTypeName*/
/*jslint sub:true*/ /* The symbols 'fromWireType' and 'toWireType' must be accessed via array notation to be closure-safe since craftInvokerFunction crafts functions as strings that can't be closured. */
var InternalError = Module.InternalError = extendError(Error, 'InternalError');
var BindingError = Module.BindingError = extendError(Error, 'BindingError');
var UnboundTypeError = Module.UnboundTypeError = extendError(BindingError, 'UnboundTypeError');
function throwInternalError(message) {
    throw new InternalError(message);
}
function throwBindingError(message) {
    throw new BindingError(message);
}
function throwUnboundTypeError(message, types) {
    var unboundTypes = [];
    var seen = {};
    function visit(type) {
        if (seen[type]) {
            return;
        }
        if (registeredTypes[type]) {
            return;
        }
        if (typeDependencies[type]) {
            typeDependencies[type].forEach(visit);
            return;
        }
        unboundTypes.push(type);
        seen[type] = true;
    }
    types.forEach(visit);
    throw new UnboundTypeError(message + ': ' + unboundTypes.map(getTypeName).join([', ']));
}
// Creates a function overload resolution table to the given method 'methodName' in the given prototype,
// if the overload table doesn't yet exist.
function ensureOverloadTable(proto, methodName, humanName) {
    if (undefined === proto[methodName].overloadTable) {
        var prevFunc = proto[methodName];
        // Inject an overload resolver function that routes to the appropriate overload based on the number of arguments.
        proto[methodName] = function() {
            // TODO This check can be removed in -O3 level "unsafe" optimizations.
            if (!proto[methodName].overloadTable.hasOwnProperty(arguments.length)) {
                throwBindingError("Function '" + humanName + "' called with an invalid number of arguments (" + arguments.length + ") - expects one of (" + proto[methodName].overloadTable + ")!");
            }
            return proto[methodName].overloadTable[arguments.length].apply(this, arguments);
        };
        // Move the previous function into the overload table.
        proto[methodName].overloadTable = [];
        proto[methodName].overloadTable[prevFunc.argCount] = prevFunc;
    }            
}
/* Registers a symbol (function, class, enum, ...) as part of the Module JS object so that
   hand-written code is able to access that symbol via 'Module.name'.
   name: The name of the symbol that's being exposed.
   value: The object itself to expose (function, class, ...)
   numArguments: For functions, specifies the number of arguments the function takes in. For other types, unused and undefined.
   To implement support for multiple overloads of a function, an 'overload selector' function is used. That selector function chooses
   the appropriate overload to call from an function overload table. This selector function is only used if multiple overloads are
   actually registered, since it carries a slight performance penalty. */
function exposePublicSymbol(name, value, numArguments) {
    if (Module.hasOwnProperty(name)) {
        if (undefined === numArguments || (undefined !== Module[name].overloadTable && undefined !== Module[name].overloadTable[numArguments])) {
            throwBindingError("Cannot register public name '" + name + "' twice");
        }
        // We are exposing a function with the same name as an existing function. Create an overload table and a function selector
        // that routes between the two.
        ensureOverloadTable(Module, name, name);
        if (Module.hasOwnProperty(numArguments)) {
            throwBindingError("Cannot register multiple overloads of a function with the same number of arguments (" + numArguments + ")!");
        }
        // Add the new function into the overload table.
        Module[name].overloadTable[numArguments] = value;
    }
    else {
        Module[name] = value;
        if (undefined !== numArguments) {
            Module[name].numArguments = numArguments;
        }
    }
}
function replacePublicSymbol(name, value, numArguments) {
    if (!Module.hasOwnProperty(name)) {
        throwInternalError('Replacing nonexistant public symbol');
    }
    // If there's an overload table for this symbol, replace the symbol in the overload table instead.
    if (undefined !== Module[name].overloadTable && undefined !== numArguments) {
        Module[name].overloadTable[numArguments] = value;
    }
    else {
        Module[name] = value;
    }
}
// from https://github.com/imvu/imvujs/blob/master/src/error.js
function extendError(baseErrorType, errorName) {
    var errorClass = createNamedFunction(errorName, function(message) {
        this.name = errorName;
        this.message = message;
        var stack = (new Error(message)).stack;
        if (stack !== undefined) {
            this.stack = this.toString() + '\n' +
                stack.replace(/^Error(:[^\n]*)?\n/, '');
        }
    });
    errorClass.prototype = Object.create(baseErrorType.prototype);
    errorClass.prototype.constructor = errorClass;
    errorClass.prototype.toString = function() {
        if (this.message === undefined) {
            return this.name;
        } else {
            return this.name + ': ' + this.message;
        }
    };
    return errorClass;
}
// from https://github.com/imvu/imvujs/blob/master/src/function.js
function createNamedFunction(name, body) {
    name = makeLegalFunctionName(name);
    /*jshint evil:true*/
    return new Function(
        "body",
        "return function " + name + "() {\n" +
        "    \"use strict\";" +
        "    return body.apply(this, arguments);\n" +
        "};\n"
    )(body);
}
function _embind_repr(v) {
    var t = typeof v;
    if (t === 'object' || t === 'array' || t === 'function') {
        return v.toString();
    } else {
        return '' + v;
    }
}
// typeID -> { toWireType: ..., fromWireType: ... }
var registeredTypes = {};
// typeID -> [callback]
var awaitingDependencies = {};
// typeID -> [dependentTypes]
var typeDependencies = {};
// class typeID -> {pointerType: ..., constPointerType: ...}
var registeredPointers = {};
function registerType(rawType, registeredInstance) {
    var name = registeredInstance.name;
    if (!rawType) {
        throwBindingError('type "' + name + '" must have a positive integer typeid pointer');
    }
    if (registeredTypes.hasOwnProperty(rawType)) {
        throwBindingError("Cannot register type '" + name + "' twice");
    }
    registeredTypes[rawType] = registeredInstance;
    delete typeDependencies[rawType];
    if (awaitingDependencies.hasOwnProperty(rawType)) {
        var callbacks = awaitingDependencies[rawType];
        delete awaitingDependencies[rawType];
        callbacks.forEach(function(cb) {
            cb();
        });
    }
}
function whenDependentTypesAreResolved(myTypes, dependentTypes, getTypeConverters) {
    myTypes.forEach(function(type) {
        typeDependencies[type] = dependentTypes;
    });
    function onComplete(typeConverters) {
        var myTypeConverters = getTypeConverters(typeConverters);
        if (myTypeConverters.length !== myTypes.length) {
            throwInternalError('Mismatched type converter count');
        }
        for (var i = 0; i < myTypes.length; ++i) {
            registerType(myTypes[i], myTypeConverters[i]);
        }
    }
    var typeConverters = new Array(dependentTypes.length);
    var unregisteredTypes = [];
    var registered = 0;
    dependentTypes.forEach(function(dt, i) {
        if (registeredTypes.hasOwnProperty(dt)) {
            typeConverters[i] = registeredTypes[dt];
        } else {
            unregisteredTypes.push(dt);
            if (!awaitingDependencies.hasOwnProperty(dt)) {
                awaitingDependencies[dt] = [];
            }
            awaitingDependencies[dt].push(function() {
                typeConverters[i] = registeredTypes[dt];
                ++registered;
                if (registered === unregisteredTypes.length) {
                    onComplete(typeConverters);
                }
            });
        }
    });
    if (0 === unregisteredTypes.length) {
        onComplete(typeConverters);
    }
}
var __charCodes = (function() {
    var codes = new Array(256);
    for (var i = 0; i < 256; ++i) {
        codes[i] = String.fromCharCode(i);
    }
    return codes;
})();
function readLatin1String(ptr) {
    var ret = "";
    var c = ptr;
    while (HEAPU8[c]) {
        ret += __charCodes[HEAPU8[c++]];
    }
    return ret;
}
function getTypeName(type) {
    var ptr = ___getTypeName(type);
    var rv = readLatin1String(ptr);
    _free(ptr);
    return rv;
}
function heap32VectorToArray(count, firstElement) {
    var array = [];
    for (var i = 0; i < count; i++) {
        array.push(HEAP32[(firstElement >> 2) + i]);
    }
    return array;
}
function requireRegisteredType(rawType, humanName) {
    var impl = registeredTypes[rawType];
    if (undefined === impl) {
        throwBindingError(humanName + " has unknown type " + getTypeName(rawType));
    }
    return impl;
}
function __embind_register_void(rawType, name) {
    name = readLatin1String(name);
    registerType(rawType, {
        name: name,
        'fromWireType': function() {
            return undefined;
        },
        'toWireType': function(destructors, o) {
            // TODO: assert if anything else is given?
            return undefined;
        },
    });
}
function __embind_register_bool(rawType, name, trueValue, falseValue) {
    name = readLatin1String(name);
    registerType(rawType, {
        name: name,
        'fromWireType': function(wt) {
            // ambiguous emscripten ABI: sometimes return values are
            // true or false, and sometimes integers (0 or 1)
            return !!wt;
        },
        'toWireType': function(destructors, o) {
            return o ? trueValue : falseValue;
        },
        destructorFunction: null, // This type does not need a destructor
    });
}
// When converting a number from JS to C++ side, the valid range of the number is
// [minRange, maxRange], inclusive.
function __embind_register_integer(primitiveType, name, minRange, maxRange) {
    name = readLatin1String(name);
    if (maxRange === -1) { // LLVM doesn't have signed and unsigned 32-bit types, so u32 literals come out as 'i32 -1'. Always treat those as max u32.
        maxRange = 4294967295;
    }
    registerType(primitiveType, {
        name: name,
        minRange: minRange,
        maxRange: maxRange,
        'fromWireType': function(value) {
            return value;
        },
        'toWireType': function(destructors, value) {
            // todo: Here we have an opportunity for -O3 level "unsafe" optimizations: we could
            // avoid the following two if()s and assume value is of proper type.
            if (typeof value !== "number" && typeof value !== "boolean") {
                throw new TypeError('Cannot convert "' + _embind_repr(value) + '" to ' + this.name);
            }
            if (value < minRange || value > maxRange) {
                throw new TypeError('Passing a number "' + _embind_repr(value) + '" from JS side to C/C++ side to an argument of type "' + name + '", which is outside the valid range [' + minRange + ', ' + maxRange + ']!');
            }
            return value | 0;
        },
        destructorFunction: null, // This type does not need a destructor
    });
}
function __embind_register_float(rawType, name) {
    name = readLatin1String(name);
    registerType(rawType, {
        name: name,
        'fromWireType': function(value) {
            return value;
        },
        'toWireType': function(destructors, value) {
            // todo: Here we have an opportunity for -O3 level "unsafe" optimizations: we could
            // avoid the following if() and assume value is of proper type.
            if (typeof value !== "number" && typeof value !== "boolean") {
                throw new TypeError('Cannot convert "' + _embind_repr(value) + '" to ' + this.name);
            }
            return value;
        },
        destructorFunction: null, // This type does not need a destructor
    });
}
function __embind_register_std_string(rawType, name) {
    name = readLatin1String(name);
    registerType(rawType, {
        name: name,
        'fromWireType': function(value) {
            var length = HEAPU32[value >> 2];
            var a = new Array(length);
            for (var i = 0; i < length; ++i) {
                a[i] = String.fromCharCode(HEAPU8[value + 4 + i]);
            }
            _free(value);
            return a.join('');
        },
        'toWireType': function(destructors, value) {
            if (value instanceof ArrayBuffer) {
                value = new Uint8Array(value);
            }
            function getTAElement(ta, index) {
                return ta[index];
            }
            function getStringElement(string, index) {
                return string.charCodeAt(index);
            }
            var getElement;
            if (value instanceof Uint8Array) {
                getElement = getTAElement;
            } else if (value instanceof Int8Array) {
                getElement = getTAElement;
            } else if (typeof value === 'string') {
                getElement = getStringElement;
            } else {
                throwBindingError('Cannot pass non-string to std::string');
            }
            // assumes 4-byte alignment
            var length = value.length;
            var ptr = _malloc(4 + length);
            HEAPU32[ptr >> 2] = length;
            for (var i = 0; i < length; ++i) {
                var charCode = getElement(value, i);
                if (charCode > 255) {
                    _free(ptr);
                    throwBindingError('String has UTF-16 code units that do not fit in 8 bits');
                }
                HEAPU8[ptr + 4 + i] = charCode;
            }
            if (destructors !== null) {
                destructors.push(_free, ptr);
            }
            return ptr;
        },
        destructorFunction: function(ptr) { _free(ptr); },
    });
}
function __embind_register_std_wstring(rawType, charSize, name) {
    name = readLatin1String(name);
    var HEAP, shift;
    if (charSize === 2) {
        HEAP = HEAPU16;
        shift = 1;
    } else if (charSize === 4) {
        HEAP = HEAPU32;
        shift = 2;
    }
    registerType(rawType, {
        name: name,
        'fromWireType': function(value) {
            var length = HEAPU32[value >> 2];
            var a = new Array(length);
            var start = (value + 4) >> shift;
            for (var i = 0; i < length; ++i) {
                a[i] = String.fromCharCode(HEAP[start + i]);
            }
            _free(value);
            return a.join('');
        },
        'toWireType': function(destructors, value) {
            // assumes 4-byte alignment
            var length = value.length;
            var ptr = _malloc(4 + length * charSize);
            HEAPU32[ptr >> 2] = length;
            var start = (ptr + 4) >> shift;
            for (var i = 0; i < length; ++i) {
                HEAP[start + i] = value.charCodeAt(i);
            }
            if (destructors !== null) {
                destructors.push(_free, ptr);
            }
            return ptr;
        },
        destructorFunction: function(ptr) { _free(ptr); },
    });
}
function __embind_register_emval(rawType, name) {
    name = readLatin1String(name);
    registerType(rawType, {
        name: name,
        'fromWireType': function(handle) {
            var rv = _emval_handle_array[handle].value;
            __emval_decref(handle);
            return rv;
        },
        'toWireType': function(destructors, value) {
            return __emval_register(value);
        },
        destructorFunction: null, // This type does not need a destructor
    });
}
function __embind_register_memory_view(rawType, name) {
    var typeMapping = [
        Int8Array,
        Uint8Array,
        Int16Array,
        Uint16Array,
        Int32Array,
        Uint32Array,
        Float32Array,
        Float64Array,        
    ];
    name = readLatin1String(name);
    registerType(rawType, {
        name: name,
        'fromWireType': function(handle) {
            var type = HEAPU32[handle >> 2];
            var size = HEAPU32[(handle >> 2) + 1]; // in elements
            var data = HEAPU32[(handle >> 2) + 2]; // byte offset into emscripten heap
            var TA = typeMapping[type];
            return new TA(HEAP8.buffer, data, size);
        },
    });
}
function runDestructors(destructors) {
    while (destructors.length) {
        var ptr = destructors.pop();
        var del = destructors.pop();
        del(ptr);
    }
}
// Function implementation of operator new, per
// http://www.ecma-international.org/publications/files/ECMA-ST/Ecma-262.pdf
// 13.2.2
// ES3
function new_(constructor, argumentList) {
    if (!(constructor instanceof Function)) {
        throw new TypeError('new_ called with constructor type ' + typeof(constructor) + " which is not a function");
    }
    /*
     * Previously, the following line was just:
     function dummy() {};
     * Unfortunately, Chrome was preserving 'dummy' as the object's name, even though at creation, the 'dummy' has the
     * correct constructor name.  Thus, objects created with IMVU.new would show up in the debugger as 'dummy', which
     * isn't very helpful.  Using IMVU.createNamedFunction addresses the issue.  Doublely-unfortunately, there's no way
     * to write a test for this behavior.  -NRD 2013.02.22
     */
    var dummy = createNamedFunction(constructor.name, function(){});
    dummy.prototype = constructor.prototype;
    var obj = new dummy;
    var r = constructor.apply(obj, argumentList);
    return (r instanceof Object) ? r : obj;
}
// The path to interop from JS code to C++ code:
// (hand-written JS code) -> (autogenerated JS invoker) -> (template-generated C++ invoker) -> (target C++ function)
// craftInvokerFunction generates the JS invoker function for each function exposed to JS through embind.
function craftInvokerFunction(humanName, argTypes, classType, cppInvokerFunc, cppTargetFunc) {
    // humanName: a human-readable string name for the function to be generated.
    // argTypes: An array that contains the embind type objects for all types in the function signature.
    //    argTypes[0] is the type object for the function return value.
    //    argTypes[1] is the type object for function this object/class type, or null if not crafting an invoker for a class method.
    //    argTypes[2...] are the actual function parameters.
    // classType: The embind type object for the class to be bound, or null if this is not a method of a class.
    // cppInvokerFunc: JS Function object to the C++-side function that interops into C++ code.
    // cppTargetFunc: Function pointer (an integer to FUNCTION_TABLE) to the target C++ function the cppInvokerFunc will end up calling.
    var argCount = argTypes.length;
    if (argCount < 2) {
        throwBindingError("argTypes array size mismatch! Must at least get return value and 'this' types!");
    }
    var isClassMethodFunc = (argTypes[1] !== null && classType !== null);
    if (!isClassMethodFunc && !FUNCTION_TABLE[cppTargetFunc]) {
        throwBindingError('Global function '+humanName+' is not defined!');
    }
    // Free functions with signature "void function()" do not need an invoker that marshalls between wire types.
// TODO: This omits argument count check - enable only at -O3 or similar.
//    if (ENABLE_UNSAFE_OPTS && argCount == 2 && argTypes[0].name == "void" && !isClassMethodFunc) {
//       return FUNCTION_TABLE[fn];
//    }
    var argsList = "";
    var argsListWired = "";
    for(var i = 0; i < argCount-2; ++i) {
        argsList += (i!==0?", ":"")+"arg"+i;
        argsListWired += (i!==0?", ":"")+"arg"+i+"Wired";
    }
    var invokerFnBody =
        "return function "+makeLegalFunctionName(humanName)+"("+argsList+") {\n" +
        "if (arguments.length !== "+(argCount - 2)+") {\n" +
            "throwBindingError('function "+humanName+" called with ' + arguments.length + ' arguments, expected "+(argCount - 2)+" args!');\n" +
        "}\n";
    // Determine if we need to use a dynamic stack to store the destructors for the function parameters.
    // TODO: Remove this completely once all function invokers are being dynamically generated.
    var needsDestructorStack = false;
    for(var i = 1; i < argTypes.length; ++i) { // Skip return value at index 0 - it's not deleted here.
        if (argTypes[i] !== null && argTypes[i].destructorFunction === undefined) { // The type does not define a destructor function - must use dynamic stack
            needsDestructorStack = true;
            break;
        }
    }
    if (needsDestructorStack) {
        invokerFnBody +=
            "var destructors = [];\n";
    }
    var dtorStack = needsDestructorStack ? "destructors" : "null";
    var args1 = ["throwBindingError", "classType", "invoker", "fn", "runDestructors", "retType", "classParam"];
    var args2 = [throwBindingError, classType, cppInvokerFunc, cppTargetFunc, runDestructors, argTypes[0], argTypes[1]];
    if (isClassMethodFunc) {
        invokerFnBody += "var thisWired = classParam.toWireType("+dtorStack+", this);\n";
    }
    for(var i = 0; i < argCount-2; ++i) {
        invokerFnBody += "var arg"+i+"Wired = argType"+i+".toWireType("+dtorStack+", arg"+i+"); // "+argTypes[i+2].name+"\n";
        args1.push("argType"+i);
        args2.push(argTypes[i+2]);
    }
    if (isClassMethodFunc) {
        argsListWired = "thisWired" + (argsListWired.length > 0 ? ", " : "") + argsListWired;
    }
    var returns = (argTypes[0].name !== "void");
    invokerFnBody +=
        (returns?"var rv = ":"") + "invoker(fn"+(argsListWired.length>0?", ":"")+argsListWired+");\n";
    if (needsDestructorStack) {
        invokerFnBody += "runDestructors(destructors);\n";
    } else {
        for(var i = isClassMethodFunc?1:2; i < argTypes.length; ++i) { // Skip return value at index 0 - it's not deleted here. Also skip class type if not a method.
            var paramName = (i === 1 ? "thisWired" : ("arg"+(i-2)+"Wired"));
            if (argTypes[i].destructorFunction !== null) {
                invokerFnBody += paramName+"_dtor("+paramName+"); // "+argTypes[i].name+"\n";
                args1.push(paramName+"_dtor");
                args2.push(argTypes[i].destructorFunction);
            }
        }
    }
    if (returns) {
        invokerFnBody += "return retType.fromWireType(rv);\n";
    }
    invokerFnBody += "}\n";
    args1.push(invokerFnBody);
    var invokerFunction = new_(Function, args1).apply(null, args2);
    return invokerFunction;
}
function __embind_register_function(name, argCount, rawArgTypesAddr, rawInvoker, fn) {
    var argTypes = heap32VectorToArray(argCount, rawArgTypesAddr);
    name = readLatin1String(name);
    rawInvoker = FUNCTION_TABLE[rawInvoker];
    exposePublicSymbol(name, function() {
        throwUnboundTypeError('Cannot call ' + name + ' due to unbound types', argTypes);
    }, argCount - 1);
    whenDependentTypesAreResolved([], argTypes, function(argTypes) {
        var invokerArgsArray = [argTypes[0] /* return value */, null /* no class 'this'*/].concat(argTypes.slice(1) /* actual params */);
        replacePublicSymbol(name, craftInvokerFunction(name, invokerArgsArray, null /* no class 'this'*/, rawInvoker, fn), argCount - 1);
        return [];
    });
}
var tupleRegistrations = {};
function __embind_register_tuple(rawType, name, rawConstructor, rawDestructor) {
    tupleRegistrations[rawType] = {
        name: readLatin1String(name),
        rawConstructor: FUNCTION_TABLE[rawConstructor],
        rawDestructor: FUNCTION_TABLE[rawDestructor],
        elements: [],
    };
}
function __embind_register_tuple_element(
    rawTupleType,
    getterReturnType,
    getter,
    getterContext,
    setterArgumentType,
    setter,
    setterContext
) {
    tupleRegistrations[rawTupleType].elements.push({
        getterReturnType: getterReturnType,
        getter: FUNCTION_TABLE[getter],
        getterContext: getterContext,
        setterArgumentType: setterArgumentType,
        setter: FUNCTION_TABLE[setter],
        setterContext: setterContext,
    });
}
function __embind_finalize_tuple(rawTupleType) {
    var reg = tupleRegistrations[rawTupleType];
    delete tupleRegistrations[rawTupleType];
    var elements = reg.elements;
    var elementsLength = elements.length;
    var elementTypes = elements.map(function(elt) { return elt.getterReturnType; }).
                concat(elements.map(function(elt) { return elt.setterArgumentType; }));
    var rawConstructor = reg.rawConstructor;
    var rawDestructor = reg.rawDestructor;
    whenDependentTypesAreResolved([rawTupleType], elementTypes, function(elementTypes) {
        elements.forEach(function(elt, i) {
            var getterReturnType = elementTypes[i];
            var getter = elt.getter;
            var getterContext = elt.getterContext;
            var setterArgumentType = elementTypes[i + elementsLength];
            var setter = elt.setter;
            var setterContext = elt.setterContext;
            elt.read = function(ptr) {
                return getterReturnType['fromWireType'](getter(getterContext, ptr));
            };
            elt.write = function(ptr, o) {
                var destructors = [];
                setter(setterContext, ptr, setterArgumentType['toWireType'](destructors, o));
                runDestructors(destructors);
            };
        });
        return [{
            name: reg.name,
            'fromWireType': function(ptr) {
                var rv = new Array(elementsLength);
                for (var i = 0; i < elementsLength; ++i) {
                    rv[i] = elements[i].read(ptr);
                }
                rawDestructor(ptr);
                return rv;
            },
            'toWireType': function(destructors, o) {
                if (elementsLength !== o.length) {
                    throw new TypeError("Incorrect number of tuple elements for " + reg.name + ": expected=" + elementsLength + ", actual=" + o.length);
                }
                var ptr = rawConstructor();
                for (var i = 0; i < elementsLength; ++i) {
                    elements[i].write(ptr, o[i]);
                }
                if (destructors !== null) {
                    destructors.push(rawDestructor, ptr);
                }
                return ptr;
            },
            destructorFunction: rawDestructor,
        }];
    });
}
var structRegistrations = {};
function __embind_register_struct(
    rawType,
    name,
    rawConstructor,
    rawDestructor
) {
    structRegistrations[rawType] = {
        name: readLatin1String(name),
        rawConstructor: FUNCTION_TABLE[rawConstructor],
        rawDestructor: FUNCTION_TABLE[rawDestructor],
        fields: [],
    };
}
function __embind_register_struct_field(
    structType,
    fieldName,
    getterReturnType,
    getter,
    getterContext,
    setterArgumentType,
    setter,
    setterContext
) {
    structRegistrations[structType].fields.push({
        fieldName: readLatin1String(fieldName),
        getterReturnType: getterReturnType,
        getter: FUNCTION_TABLE[getter],
        getterContext: getterContext,
        setterArgumentType: setterArgumentType,
        setter: FUNCTION_TABLE[setter],
        setterContext: setterContext,
    });
}
function __embind_finalize_struct(structType) {
    var reg = structRegistrations[structType];
    delete structRegistrations[structType];
    var rawConstructor = reg.rawConstructor;
    var rawDestructor = reg.rawDestructor;
    var fieldRecords = reg.fields;
    var fieldTypes = fieldRecords.map(function(field) { return field.getterReturnType; }).
              concat(fieldRecords.map(function(field) { return field.setterArgumentType; }));
    whenDependentTypesAreResolved([structType], fieldTypes, function(fieldTypes) {
        var fields = {};
        fieldRecords.forEach(function(field, i) {
            var fieldName = field.fieldName;
            var getterReturnType = fieldTypes[i];
            var getter = field.getter;
            var getterContext = field.getterContext;
            var setterArgumentType = fieldTypes[i + fieldRecords.length];
            var setter = field.setter;
            var setterContext = field.setterContext;
            fields[fieldName] = {
                read: function(ptr) {
                    return getterReturnType['fromWireType'](
                        getter(getterContext, ptr));
                },
                write: function(ptr, o) {
                    var destructors = [];
                    setter(setterContext, ptr, setterArgumentType['toWireType'](destructors, o));
                    runDestructors(destructors);
                }
            };
        });
        return [{
            name: reg.name,
            'fromWireType': function(ptr) {
                var rv = {};
                for (var i in fields) {
                    rv[i] = fields[i].read(ptr);
                }
                rawDestructor(ptr);
                return rv;
            },
            'toWireType': function(destructors, o) {
                // todo: Here we have an opportunity for -O3 level "unsafe" optimizations:
                // assume all fields are present without checking.
                for (var fieldName in fields) {
                    if (!(fieldName in o)) {
                        throw new TypeError('Missing field');
                    }
                }
                var ptr = rawConstructor();
                for (fieldName in fields) {
                    fields[fieldName].write(ptr, o[fieldName]);
                }
                if (destructors !== null) {
                    destructors.push(rawDestructor, ptr);
                }
                return ptr;
            },
            destructorFunction: rawDestructor,
        }];
    });
}
var genericPointerToWireType = function(destructors, handle) {
    if (handle === null) {
        if (this.isReference) {
            throwBindingError('null is not a valid ' + this.name);
        }
        if (this.isSmartPointer) {
            var ptr = this.rawConstructor();
            if (destructors !== null) {
                destructors.push(this.rawDestructor, ptr);
            }
            return ptr;
        } else {
            return 0;
        }
    }
    if (!handle.$$) {
        throwBindingError('Cannot pass "' + _embind_repr(handle) + '" as a ' + this.name);
    }
    if (!handle.$$.ptr) {
        throwBindingError('Cannot pass deleted object as a pointer of type ' + this.name);
    }
    if (!this.isConst && handle.$$.ptrType.isConst) {
        throwBindingError('Cannot convert argument of type ' + (handle.$$.smartPtrType ? handle.$$.smartPtrType.name : handle.$$.ptrType.name) + ' to parameter type ' + this.name);
    }
    var handleClass = handle.$$.ptrType.registeredClass;
    var ptr = upcastPointer(handle.$$.ptr, handleClass, this.registeredClass);
    if (this.isSmartPointer) {
        // TODO: this is not strictly true
        // We could support BY_EMVAL conversions from raw pointers to smart pointers
        // because the smart pointer can hold a reference to the handle
        if (undefined === handle.$$.smartPtr) {
            throwBindingError('Passing raw pointer to smart pointer is illegal');
        }
        switch (this.sharingPolicy) {
            case 0: // NONE
                // no upcasting
                if (handle.$$.smartPtrType === this) {
                    ptr = handle.$$.smartPtr;
                } else {
                    throwBindingError('Cannot convert argument of type ' + (handle.$$.smartPtrType ? handle.$$.smartPtrType.name : handle.$$.ptrType.name) + ' to parameter type ' + this.name);
                }
                break;
            case 1: // INTRUSIVE
                ptr = handle.$$.smartPtr;
                break;
            case 2: // BY_EMVAL
                if (handle.$$.smartPtrType === this) {
                    ptr = handle.$$.smartPtr;
                } else {
                    var clonedHandle = handle.clone();
                    ptr = this.rawShare(
                        ptr,
                        __emval_register(function() {
                            clonedHandle.delete();
                        })
                    );
                    if (destructors !== null) {
                        destructors.push(this.rawDestructor, ptr);
                    }
                }
                break;
            default:
                throwBindingError('Unsupporting sharing policy');
        }
    }
    return ptr;
};
// If we know a pointer type is not going to have SmartPtr logic in it, we can
// special-case optimize it a bit (compare to genericPointerToWireType)
var constNoSmartPtrRawPointerToWireType = function(destructors, handle) {
    if (handle === null) {
        if (this.isReference) {
            throwBindingError('null is not a valid ' + this.name);
        }
        return 0;
    }
    if (!handle.$$) {
        throwBindingError('Cannot pass "' + _embind_repr(handle) + '" as a ' + this.name);
    }
    if (!handle.$$.ptr) {
        throwBindingError('Cannot pass deleted object as a pointer of type ' + this.name);
    }
    var handleClass = handle.$$.ptrType.registeredClass;
    var ptr = upcastPointer(handle.$$.ptr, handleClass, this.registeredClass);
    return ptr;
};
// An optimized version for non-const method accesses - there we must additionally restrict that
// the pointer is not a const-pointer.
var nonConstNoSmartPtrRawPointerToWireType = function(destructors, handle) {
    if (handle === null) {
        if (this.isReference) {
            throwBindingError('null is not a valid ' + this.name);
        }
        return 0;
    }
    if (!handle.$$) {
        throwBindingError('Cannot pass "' + _embind_repr(handle) + '" as a ' + this.name);
    }
    if (!handle.$$.ptr) {
        throwBindingError('Cannot pass deleted object as a pointer of type ' + this.name);
    }
    if (handle.$$.ptrType.isConst) {
        throwBindingError('Cannot convert argument of type ' + handle.$$.ptrType.name + ' to parameter type ' + this.name);
    }
    var handleClass = handle.$$.ptrType.registeredClass;
    var ptr = upcastPointer(handle.$$.ptr, handleClass, this.registeredClass);
    return ptr;
};
function RegisteredPointer(
    name,
    registeredClass,
    isReference,
    isConst,
    // smart pointer properties
    isSmartPointer,
    pointeeType,
    sharingPolicy,
    rawGetPointee,
    rawConstructor,
    rawShare,
    rawDestructor
) {
    this.name = name;
    this.registeredClass = registeredClass;
    this.isReference = isReference;
    this.isConst = isConst;
    // smart pointer properties
    this.isSmartPointer = isSmartPointer;
    this.pointeeType = pointeeType;
    this.sharingPolicy = sharingPolicy;
    this.rawGetPointee = rawGetPointee;
    this.rawConstructor = rawConstructor;
    this.rawShare = rawShare;
    this.rawDestructor = rawDestructor;
    if (!isSmartPointer && registeredClass.baseClass === undefined) {
        if (isConst) {
            this['toWireType'] = constNoSmartPtrRawPointerToWireType;
            this.destructorFunction = null;
        } else {
            this['toWireType'] = nonConstNoSmartPtrRawPointerToWireType;
            this.destructorFunction = null;
        }
    } else {
        this['toWireType'] = genericPointerToWireType;
        // Here we must leave this.destructorFunction undefined, since whether genericPointerToWireType returns
        // a pointer that needs to be freed up is runtime-dependent, and cannot be evaluated at registration time.
        // TODO: Create an alternative mechanism that allows removing the use of var destructors = []; array in 
        //       craftInvokerFunction altogether.
    }
}
RegisteredPointer.prototype.getPointee = function(ptr) {
    if (this.rawGetPointee) {
        ptr = this.rawGetPointee(ptr);
    }
    return ptr;
};
RegisteredPointer.prototype.destructor = function(ptr) {
    if (this.rawDestructor) {
        this.rawDestructor(ptr);
    }
};
RegisteredPointer.prototype['fromWireType'] = function(ptr) {
    // ptr is a raw pointer (or a raw smartpointer)
    // rawPointer is a maybe-null raw pointer
    var rawPointer = this.getPointee(ptr);
    if (!rawPointer) {
        this.destructor(ptr);
        return null;
    }
    function makeDefaultHandle() {
        if (this.isSmartPointer) {
            return makeClassHandle(this.registeredClass.instancePrototype, {
                ptrType: this.pointeeType,
                ptr: rawPointer,
                smartPtrType: this,
                smartPtr: ptr,
            });
        } else {
            return makeClassHandle(this.registeredClass.instancePrototype, {
                ptrType: this,
                ptr: ptr,
            });
        }
    }
    var actualType = this.registeredClass.getActualType(rawPointer);
    var registeredPointerRecord = registeredPointers[actualType];
    if (!registeredPointerRecord) {
        return makeDefaultHandle.call(this);
    }
    var toType;
    if (this.isConst) {
        toType = registeredPointerRecord.constPointerType;
    } else {
        toType = registeredPointerRecord.pointerType;
    }
    var dp = downcastPointer(
        rawPointer,
        this.registeredClass,
        toType.registeredClass);
    if (dp === null) {
        return makeDefaultHandle.call(this);
    }
    if (this.isSmartPointer) {
        return makeClassHandle(toType.registeredClass.instancePrototype, {
            ptrType: toType,
            ptr: dp,
            smartPtrType: this,
            smartPtr: ptr,
        });
    } else {
        return makeClassHandle(toType.registeredClass.instancePrototype, {
            ptrType: toType,
            ptr: dp,
        });
    }
};
function makeClassHandle(prototype, record) {
    if (!record.ptrType || !record.ptr) {
        throwInternalError('makeClassHandle requires ptr and ptrType');
    }
    var hasSmartPtrType = !!record.smartPtrType;
    var hasSmartPtr = !!record.smartPtr;
    if (hasSmartPtrType !== hasSmartPtr) {
        throwInternalError('Both smartPtrType and smartPtr must be specified');
    }
    record.count = { value: 1 };
    return Object.create(prototype, {
        $$: {
            value: record,
        },
    });
}
// root of all pointer and smart pointer handles in embind
function ClassHandle() {
}
function getInstanceTypeName(handle) {
    return handle.$$.ptrType.registeredClass.name;
}
ClassHandle.prototype.isAliasOf = function(other) {
    if (!(this instanceof ClassHandle)) {
        return false;
    }
    if (!(other instanceof ClassHandle)) {
        return false;
    }
    var leftClass = this.$$.ptrType.registeredClass;
    var left = this.$$.ptr;
    var rightClass = other.$$.ptrType.registeredClass;
    var right = other.$$.ptr;
    while (leftClass.baseClass) {
        left = leftClass.upcast(left);
        leftClass = leftClass.baseClass;
    }
    while (rightClass.baseClass) {
        right = rightClass.upcast(right);
        rightClass = rightClass.baseClass;
    }
    return leftClass === rightClass && left === right;
};
function throwInstanceAlreadyDeleted(obj) {
    throwBindingError(getInstanceTypeName(obj) + ' instance already deleted');
}
ClassHandle.prototype.clone = function() {
    if (!this.$$.ptr) {
        throwInstanceAlreadyDeleted(this);
    }
    var clone = Object.create(Object.getPrototypeOf(this), {
        $$: {
            value: shallowCopy(this.$$),
        }
    });
    clone.$$.count.value += 1;
    return clone;
};
function runDestructor(handle) {
    var $$ = handle.$$;
    if ($$.smartPtr) {
        $$.smartPtrType.rawDestructor($$.smartPtr);
    } else {
        $$.ptrType.registeredClass.rawDestructor($$.ptr);
    }
}
ClassHandle.prototype['delete'] = function ClassHandle_delete() {
    if (!this.$$.ptr) {
        throwInstanceAlreadyDeleted(this);
    }
    if (this.$$.deleteScheduled) {
        throwBindingError('Object already scheduled for deletion');
    }
    this.$$.count.value -= 1;
    if (0 === this.$$.count.value) {
        runDestructor(this);
    }
    this.$$.smartPtr = undefined;
    this.$$.ptr = undefined;
};
var deletionQueue = [];
ClassHandle.prototype['isDeleted'] = function isDeleted() {
    return !this.$$.ptr;
};
ClassHandle.prototype['deleteLater'] = function deleteLater() {
    if (!this.$$.ptr) {
        throwInstanceAlreadyDeleted(this);
    }
    if (this.$$.deleteScheduled) {
        throwBindingError('Object already scheduled for deletion');
    }
    deletionQueue.push(this);
    if (deletionQueue.length === 1 && delayFunction) {
        delayFunction(flushPendingDeletes);
    }
    this.$$.deleteScheduled = true;
    return this;
};
function flushPendingDeletes() {
    while (deletionQueue.length) {
        var obj = deletionQueue.pop();
        obj.$$.deleteScheduled = false;
        obj['delete']();
    }
}
Module['flushPendingDeletes'] = flushPendingDeletes;
var delayFunction;
Module['setDelayFunction'] = function setDelayFunction(fn) {
    delayFunction = fn;
    if (deletionQueue.length && delayFunction) {
        delayFunction(flushPendingDeletes);
    }
};
function RegisteredClass(
    name,
    constructor,
    instancePrototype,
    rawDestructor,
    baseClass,
    getActualType,
    upcast,
    downcast
) {
    this.name = name;
    this.constructor = constructor;
    this.instancePrototype = instancePrototype;
    this.rawDestructor = rawDestructor;
    this.baseClass = baseClass;
    this.getActualType = getActualType;
    this.upcast = upcast;
    this.downcast = downcast;
}
function shallowCopy(o) {
    var rv = {};
    for (var k in o) {
        rv[k] = o[k];
    }
    return rv;
}
function __embind_register_class(
    rawType,
    rawPointerType,
    rawConstPointerType,
    baseClassRawType,
    getActualType,
    upcast,
    downcast,
    name,
    rawDestructor
) {
    name = readLatin1String(name);
    rawDestructor = FUNCTION_TABLE[rawDestructor];
    getActualType = FUNCTION_TABLE[getActualType];
    upcast = FUNCTION_TABLE[upcast];
    downcast = FUNCTION_TABLE[downcast];
    var legalFunctionName = makeLegalFunctionName(name);
    exposePublicSymbol(legalFunctionName, function() {
        // this code cannot run if baseClassRawType is zero
        throwUnboundTypeError('Cannot construct ' + name + ' due to unbound types', [baseClassRawType]);
    });
    whenDependentTypesAreResolved(
        [rawType, rawPointerType, rawConstPointerType],
        baseClassRawType ? [baseClassRawType] : [],
        function(base) {
            base = base[0];
            var baseClass;
            var basePrototype;
            if (baseClassRawType) {
                baseClass = base.registeredClass;
                basePrototype = baseClass.instancePrototype;
            } else {
                basePrototype = ClassHandle.prototype;
            }
            var constructor = createNamedFunction(legalFunctionName, function() {
                if (Object.getPrototypeOf(this) !== instancePrototype) {
                    throw new BindingError("Use 'new' to construct " + name);
                }
                if (undefined === registeredClass.constructor_body) {
                    throw new BindingError(name + " has no accessible constructor");
                }
                var body = registeredClass.constructor_body[arguments.length];
                if (undefined === body) {
                    throw new BindingError("Tried to invoke ctor of " + name + " with invalid number of parameters (" + arguments.length + ") - expected (" + Object.keys(registeredClass.constructor_body).toString() + ") parameters instead!");
                }
                return body.apply(this, arguments);
            });
            var instancePrototype = Object.create(basePrototype, {
                constructor: { value: constructor },
            });
            constructor.prototype = instancePrototype;
            var registeredClass = new RegisteredClass(
                name,
                constructor,
                instancePrototype,
                rawDestructor,
                baseClass,
                getActualType,
                upcast,
                downcast);
            var referenceConverter = new RegisteredPointer(
                name,
                registeredClass,
                true,
                false,
                false);
            var pointerConverter = new RegisteredPointer(
                name + '*',
                registeredClass,
                false,
                false,
                false);
            var constPointerConverter = new RegisteredPointer(
                name + ' const*',
                registeredClass,
                false,
                true,
                false);
            registeredPointers[rawType] = {
                pointerType: pointerConverter,
                constPointerType: constPointerConverter
            };
            replacePublicSymbol(legalFunctionName, constructor);
            return [referenceConverter, pointerConverter, constPointerConverter];
        }
    );
}
function __embind_register_class_constructor(
    rawClassType,
    argCount,
    rawArgTypesAddr,
    invoker,
    rawConstructor
) {
    var rawArgTypes = heap32VectorToArray(argCount, rawArgTypesAddr);
    invoker = FUNCTION_TABLE[invoker];
    whenDependentTypesAreResolved([], [rawClassType], function(classType) {
        classType = classType[0];
        var humanName = 'constructor ' + classType.name;
        if (undefined === classType.registeredClass.constructor_body) {
            classType.registeredClass.constructor_body = [];
        }
        if (undefined !== classType.registeredClass.constructor_body[argCount - 1]) {
            throw new BindingError("Cannot register multiple constructors with identical number of parameters (" + (argCount-1) + ") for class '" + classType.name + "'! Overload resolution is currently only performed using the parameter count, not actual type info!");
        }
        classType.registeredClass.constructor_body[argCount - 1] = function() {
            throwUnboundTypeError('Cannot construct ' + classType.name + ' due to unbound types', rawArgTypes);
        };
        whenDependentTypesAreResolved([], rawArgTypes, function(argTypes) {
            classType.registeredClass.constructor_body[argCount - 1] = function() {
                if (arguments.length !== argCount - 1) {
                    throwBindingError(humanName + ' called with ' + arguments.length + ' arguments, expected ' + (argCount-1));
                }
                var destructors = [];
                var args = new Array(argCount);
                args[0] = rawConstructor;
                for (var i = 1; i < argCount; ++i) {
                    args[i] = argTypes[i]['toWireType'](destructors, arguments[i - 1]);
                }
                var ptr = invoker.apply(null, args);
                runDestructors(destructors);
                return argTypes[0]['fromWireType'](ptr);
            };
            return [];
        });
        return [];
    });
}
function downcastPointer(ptr, ptrClass, desiredClass) {
    if (ptrClass === desiredClass) {
        return ptr;
    }
    if (undefined === desiredClass.baseClass) {
        return null; // no conversion
    }
    // O(depth) stack space used
    return desiredClass.downcast(
        downcastPointer(ptr, ptrClass, desiredClass.baseClass));
}
function upcastPointer(ptr, ptrClass, desiredClass) {
    while (ptrClass !== desiredClass) {
        if (!ptrClass.upcast) {
            throwBindingError("Expected null or instance of " + desiredClass.name + ", got an instance of " + ptrClass.name);
        }
        ptr = ptrClass.upcast(ptr);
        ptrClass = ptrClass.baseClass;
    }
    return ptr;
}
function validateThis(this_, classType, humanName) {
    if (!(this_ instanceof Object)) {
        throwBindingError(humanName + ' with invalid "this": ' + this_);
    }
    if (!(this_ instanceof classType.registeredClass.constructor)) {
        throwBindingError(humanName + ' incompatible with "this" of type ' + this_.constructor.name);
    }
    if (!this_.$$.ptr) {
        throwBindingError('cannot call emscripten binding method ' + humanName + ' on deleted object');
    }
    // todo: kill this
    return upcastPointer(
        this_.$$.ptr,
        this_.$$.ptrType.registeredClass,
        classType.registeredClass);
}
function __embind_register_class_function(
    rawClassType,
    methodName,
    argCount,
    rawArgTypesAddr, // [ReturnType, ThisType, Args...]
    rawInvoker,
    context
) {
    var rawArgTypes = heap32VectorToArray(argCount, rawArgTypesAddr);
    methodName = readLatin1String(methodName);
    rawInvoker = FUNCTION_TABLE[rawInvoker];
    whenDependentTypesAreResolved([], [rawClassType], function(classType) {
        classType = classType[0];
        var humanName = classType.name + '.' + methodName;
        var unboundTypesHandler = function() {
            throwUnboundTypeError('Cannot call ' + humanName + ' due to unbound types', rawArgTypes);
        };
        var proto = classType.registeredClass.instancePrototype;
        var method = proto[methodName];
        if (undefined === method || (undefined === method.overloadTable && method.className !== classType.name && method.argCount === argCount-2)) {
            // This is the first overload to be registered, OR we are replacing a function in the base class with a function in the derived class.
            unboundTypesHandler.argCount = argCount-2;
            unboundTypesHandler.className = classType.name;
            proto[methodName] = unboundTypesHandler;
        } else {
            // There was an existing function with the same name registered. Set up a function overload routing table.
            ensureOverloadTable(proto, methodName, humanName);
            proto[methodName].overloadTable[argCount-2] = unboundTypesHandler;
        }
        whenDependentTypesAreResolved([], rawArgTypes, function(argTypes) {
            var memberFunction = craftInvokerFunction(humanName, argTypes, classType, rawInvoker, context);
            // Replace the initial unbound-handler-stub function with the appropriate member function, now that all types
            // are resolved. If multiple overloads are registered for this function, the function goes into an overload table.
            if (undefined === proto[methodName].overloadTable) {
                proto[methodName] = memberFunction;
            } else {
                proto[methodName].overloadTable[argCount-2] = memberFunction;
            }
            return [];
        });
        return [];
    });
}
function __embind_register_class_class_function(
    rawClassType,
    methodName,
    argCount,
    rawArgTypesAddr,
    rawInvoker,
    fn
) {
    var rawArgTypes = heap32VectorToArray(argCount, rawArgTypesAddr);
    methodName = readLatin1String(methodName);
    rawInvoker = FUNCTION_TABLE[rawInvoker];
    whenDependentTypesAreResolved([], [rawClassType], function(classType) {
        classType = classType[0];
        var humanName = classType.name + '.' + methodName;
        var unboundTypesHandler = function() {
                throwUnboundTypeError('Cannot call ' + humanName + ' due to unbound types', rawArgTypes);
            };
        var proto = classType.registeredClass.constructor;
        if (undefined === proto[methodName]) {
            // This is the first function to be registered with this name.
            unboundTypesHandler.argCount = argCount-1;
            proto[methodName] = unboundTypesHandler;
        } else {
            // There was an existing function with the same name registered. Set up a function overload routing table.
            ensureOverloadTable(proto, methodName, humanName);
            proto[methodName].overloadTable[argCount-1] = unboundTypesHandler;
        }
        whenDependentTypesAreResolved([], rawArgTypes, function(argTypes) {
            // Replace the initial unbound-types-handler stub with the proper function. If multiple overloads are registered,
            // the function handlers go into an overload table.
            var invokerArgsArray = [argTypes[0] /* return value */, null /* no class 'this'*/].concat(argTypes.slice(1) /* actual params */);
            var func = craftInvokerFunction(humanName, invokerArgsArray, null /* no class 'this'*/, rawInvoker, fn);
            if (undefined === proto[methodName].overloadTable) {
                proto[methodName] = func;
            } else {
                proto[methodName].overloadTable[argCount-1] = func;
            }
            return [];
        });
        return [];
    });
}
function __embind_register_class_property(
    classType,
    fieldName,
    getterReturnType,
    getter,
    getterContext,
    setterArgumentType,
    setter,
    setterContext
) {
    fieldName = readLatin1String(fieldName);
    getter = FUNCTION_TABLE[getter];
    whenDependentTypesAreResolved([], [classType], function(classType) {
        classType = classType[0];
        var humanName = classType.name + '.' + fieldName;
        var desc = {
            get: function() {
                throwUnboundTypeError('Cannot access ' + humanName + ' due to unbound types', [getterReturnType, setterArgumentType]);
            },
            enumerable: true,
            configurable: true
        };
        if (setter) {
            desc.set = function() {
                throwUnboundTypeError('Cannot access ' + humanName + ' due to unbound types', [getterReturnType, setterArgumentType]);
            };
        } else {
            desc.set = function(v) {
                throwBindingError(humanName + ' is a read-only property');
            };
        }
        Object.defineProperty(classType.registeredClass.instancePrototype, fieldName, desc);
        whenDependentTypesAreResolved(
            [],
            (setter ? [getterReturnType, setterArgumentType] : [getterReturnType]),
        function(types) {
            var getterReturnType = types[0];
            var desc = {
                get: function() {
                    var ptr = validateThis(this, classType, humanName + ' getter');
                    return getterReturnType['fromWireType'](getter(getterContext, ptr));
                },
                enumerable: true
            };
            if (setter) {
                setter = FUNCTION_TABLE[setter];
                var setterArgumentType = types[1];
                desc.set = function(v) {
                    var ptr = validateThis(this, classType, humanName + ' setter');
                    var destructors = [];
                    setter(setterContext, ptr, setterArgumentType['toWireType'](destructors, v));
                    runDestructors(destructors);
                };
            }
            Object.defineProperty(classType.registeredClass.instancePrototype, fieldName, desc);
            return [];
        });
        return [];
    });
}
var char_0 = '0'.charCodeAt(0);
var char_9 = '9'.charCodeAt(0);
function makeLegalFunctionName(name) {
    name = name.replace(/[^a-zA-Z0-9_]/g, '$');
    var f = name.charCodeAt(0);
    if (f >= char_0 && f <= char_9) {
        return '_' + name;
    } else {
        return name;
    }
}
function __embind_register_smart_ptr(
    rawType,
    rawPointeeType,
    name,
    sharingPolicy,
    rawGetPointee,
    rawConstructor,
    rawShare,
    rawDestructor
) {
    name = readLatin1String(name);
    rawGetPointee = FUNCTION_TABLE[rawGetPointee];
    rawConstructor = FUNCTION_TABLE[rawConstructor];
    rawShare = FUNCTION_TABLE[rawShare];
    rawDestructor = FUNCTION_TABLE[rawDestructor];
    whenDependentTypesAreResolved([rawType], [rawPointeeType], function(pointeeType) {
        pointeeType = pointeeType[0];
        var registeredPointer = new RegisteredPointer(
            name,
            pointeeType.registeredClass,
            false,
            false,
            // smart pointer properties
            true,
            pointeeType,
            sharingPolicy,
            rawGetPointee,
            rawConstructor,
            rawShare,
            rawDestructor);
        return [registeredPointer];
    });
}
function __embind_register_enum(
    rawType,
    name
) {
    name = readLatin1String(name);
    function constructor() {
    }
    constructor.values = {};
    registerType(rawType, {
        name: name,
        constructor: constructor,
        'fromWireType': function(c) {
            return this.constructor.values[c];
        },
        'toWireType': function(destructors, c) {
            return c.value;
        },
        destructorFunction: null,
    });
    exposePublicSymbol(name, constructor);
}
function __embind_register_enum_value(
    rawEnumType,
    name,
    enumValue
) {
    var enumType = requireRegisteredType(rawEnumType, 'enum');
    name = readLatin1String(name);
    var Enum = enumType.constructor;
    var Value = Object.create(enumType.constructor.prototype, {
        value: {value: enumValue},
        constructor: {value: createNamedFunction(enumType.name + '_' + name, function() {})},
    });
    Enum.values[enumValue] = Value;
    Enum[name] = Value;
}
function __embind_register_constant(name, type, value) {
    name = readLatin1String(name);
    whenDependentTypesAreResolved([], [type], function(type) {
        type = type[0];
        Module[name] = type['fromWireType'](value);
        return [];
    });
}
/*global Module:true, Runtime*/
/*global HEAP32*/
/*global new_*/
/*global createNamedFunction*/
/*global readLatin1String, writeStringToMemory*/
/*global requireRegisteredType, throwBindingError*/
var Module = Module || {};
var _emval_handle_array = [{}]; // reserve zero
var _emval_free_list = [];
// Public JS API
/** @expose */
Module.count_emval_handles = function() {
    var count = 0;
    for (var i = 1; i < _emval_handle_array.length; ++i) {
        if (_emval_handle_array[i] !== undefined) {
            ++count;
        }
    }
    return count;
};
/** @expose */
Module.get_first_emval = function() {
    for (var i = 1; i < _emval_handle_array.length; ++i) {
        if (_emval_handle_array[i] !== undefined) {
            return _emval_handle_array[i];
        }
    }
    return null;
};
// Private C++ API
var _emval_symbols = {}; // address -> string
function __emval_register_symbol(address) {
    _emval_symbols[address] = readLatin1String(address);
}
function getStringOrSymbol(address) {
    var symbol = _emval_symbols[address];
    if (symbol === undefined) {
        return readLatin1String(address);
    } else {
        return symbol;
    }
}
function requireHandle(handle) {
    if (!handle) {
        throwBindingError('Cannot use deleted val. handle = ' + handle);
    }
}
function __emval_register(value) {
    var handle = _emval_free_list.length ?
        _emval_free_list.pop() :
        _emval_handle_array.length;
    _emval_handle_array[handle] = {refcount: 1, value: value};
    return handle;
}
function __emval_incref(handle) {
    if (handle) {
        _emval_handle_array[handle].refcount += 1;
    }
}
function __emval_decref(handle) {
    if (handle && 0 === --_emval_handle_array[handle].refcount) {
        _emval_handle_array[handle] = undefined;
        _emval_free_list.push(handle);
    }
}
function __emval_new_array() {
    return __emval_register([]);
}
function __emval_new_object() {
    return __emval_register({});
}
function __emval_undefined() {
    return __emval_register(undefined);
}
function __emval_null() {
    return __emval_register(null);
}
function __emval_new_cstring(v) {
    return __emval_register(getStringOrSymbol(v));
}
function __emval_take_value(type, v) {
    type = requireRegisteredType(type, '_emval_take_value');
    v = type.fromWireType(v);
    return __emval_register(v);
}
var __newers = {}; // arity -> function
function craftEmvalAllocator(argCount) {
    /*This function returns a new function that looks like this:
    function emval_allocator_3(handle, argTypes, arg0Wired, arg1Wired, arg2Wired) {
        var argType0 = requireRegisteredType(HEAP32[(argTypes >> 2)], "parameter 0");
        var arg0 = argType0.fromWireType(arg0Wired);
        var argType1 = requireRegisteredType(HEAP32[(argTypes >> 2) + 1], "parameter 1");
        var arg1 = argType1.fromWireType(arg1Wired);
        var argType2 = requireRegisteredType(HEAP32[(argTypes >> 2) + 2], "parameter 2");
        var arg2 = argType2.fromWireType(arg2Wired);
        var constructor = _emval_handle_array[handle].value;
        var emval = new constructor(arg0, arg1, arg2);
        return emval;
    } */
    var args1 = ["requireRegisteredType", "HEAP32", "_emval_handle_array", "__emval_register"];
    var args2 = [requireRegisteredType, HEAP32, _emval_handle_array, __emval_register];
    var argsList = "";
    var argsListWired = "";
    for(var i = 0; i < argCount; ++i) {
        argsList += (i!==0?", ":"")+"arg"+i; // 'arg0, arg1, ..., argn'
        argsListWired += ", arg"+i+"Wired"; // ', arg0Wired, arg1Wired, ..., argnWired'
    }
    var invokerFnBody =
        "return function emval_allocator_"+argCount+"(handle, argTypes " + argsListWired + ") {\n";
    for(var i = 0; i < argCount; ++i) {
        invokerFnBody += 
            "var argType"+i+" = requireRegisteredType(HEAP32[(argTypes >> 2) + "+i+"], \"parameter "+i+"\");\n" +
            "var arg"+i+" = argType"+i+".fromWireType(arg"+i+"Wired);\n";
    }
    invokerFnBody +=
        "var constructor = _emval_handle_array[handle].value;\n" +
        "var obj = new constructor("+argsList+");\n" +
        "return __emval_register(obj);\n" +
        "}\n";
    args1.push(invokerFnBody);
    var invokerFunction = new_(Function, args1).apply(null, args2);
    return invokerFunction;
}
function __emval_new(handle, argCount, argTypes) {
    requireHandle(handle);
    var newer = __newers[argCount];
    if (!newer) {
        newer = craftEmvalAllocator(argCount);
        __newers[argCount] = newer;
    }
    if (argCount === 0) {
        return newer(handle, argTypes);
    } else if (argCount === 1) {
        return newer(handle, argTypes, arguments[3]);
    } else if (argCount === 2) {
        return newer(handle, argTypes, arguments[3], arguments[4]);
    } else if (argCount === 3) {
        return newer(handle, argTypes, arguments[3], arguments[4], arguments[5]);
    } else if (argCount === 4) {
        return newer(handle, argTypes, arguments[3], arguments[4], arguments[5], arguments[6]);
    } else {
        // This is a slow path! (.apply and .splice are slow), so a few specializations are present above.
        return newer.apply(null, arguments.splice(1));
    }
}
// appease jshint (technically this code uses eval)
var global = (function(){return Function;})()('return this')();
function __emval_get_global(name) {
    name = getStringOrSymbol(name);
    return __emval_register(global[name]);
}
function __emval_get_module_property(name) {
    name = getStringOrSymbol(name);
    return __emval_register(Module[name]);
}
function __emval_get_property(handle, key) {
    requireHandle(handle);
    return __emval_register(_emval_handle_array[handle].value[_emval_handle_array[key].value]);
}
function __emval_set_property(handle, key, value) {
    requireHandle(handle);
    _emval_handle_array[handle].value[_emval_handle_array[key].value] = _emval_handle_array[value].value;
}
function __emval_as(handle, returnType) {
    requireHandle(handle);
    returnType = requireRegisteredType(returnType, 'emval::as');
    var destructors = [];
    // caller owns destructing
    return returnType.toWireType(destructors, _emval_handle_array[handle].value);
}
function parseParameters(argCount, argTypes, argWireTypes) {
    var a = new Array(argCount);
    for (var i = 0; i < argCount; ++i) {
        var argType = requireRegisteredType(
            HEAP32[(argTypes >> 2) + i],
            "parameter " + i);
        a[i] = argType.fromWireType(argWireTypes[i]);
    }
    return a;
}
function __emval_call(handle, argCount, argTypes) {
    requireHandle(handle);
    var types = lookupTypes(argCount, argTypes);
    var args = new Array(argCount);
    for (var i = 0; i < argCount; ++i) {
        args[i] = types[i].fromWireType(arguments[3 + i]);
    }
    var fn = _emval_handle_array[handle].value;
    var rv = fn.apply(undefined, args);
    return __emval_register(rv);
}
function lookupTypes(argCount, argTypes, argWireTypes) {
    var a = new Array(argCount);
    for (var i = 0; i < argCount; ++i) {
        a[i] = requireRegisteredType(
            HEAP32[(argTypes >> 2) + i],
            "parameter " + i);
    }
    return a;
}
function __emval_get_method_caller(argCount, argTypes) {
    var types = lookupTypes(argCount, argTypes);
    var retType = types[0];
    var signatureName = retType.name + "_$" + types.slice(1).map(function (t) { return t.name; }).join("_") + "$";
    var args1 = ["Runtime", "createNamedFunction", "requireHandle", "getStringOrSymbol", "_emval_handle_array", "retType"];
    var args2 = [Runtime, createNamedFunction, requireHandle, getStringOrSymbol, _emval_handle_array, retType];
    var argsList = ""; // 'arg0, arg1, arg2, ... , argN'
    var argsListWired = ""; // 'arg0Wired, ..., argNWired'
    for (var i = 0; i < argCount - 1; ++i) {
        argsList += (i !== 0 ? ", " : "") + "arg" + i;
        argsListWired += ", arg" + i + "Wired";
        args1.push("argType" + i);
        args2.push(types[1 + i]);
    }
    var invokerFnBody =
        "return Runtime.addFunction(createNamedFunction('" + signatureName + "', function (handle, name" + argsListWired + ") {\n" +
        "requireHandle(handle);\n" +
        "name = getStringOrSymbol(name);\n";
    for (var i = 0; i < argCount - 1; ++i) {
        invokerFnBody += "var arg" + i + " = argType" + i + ".fromWireType(arg" + i + "Wired);\n";
    }
    invokerFnBody +=
        "var obj = _emval_handle_array[handle].value;\n" +
        "return retType.toWireType(null, obj[name](" + argsList + "));\n" + 
        "}));\n";
    args1.push(invokerFnBody);
    var invokerFunction = new_(Function, args1).apply(null, args2);
    return invokerFunction;
}
function __emval_has_function(handle, name) {
    name = getStringOrSymbol(name);
    return _emval_handle_array[handle].value[name] instanceof Function;
}
if (Module['preInit']) {
  if (typeof Module['preInit'] == 'function') Module['preInit'] = [Module['preInit']];
  while (Module['preInit'].length > 0) {
    Module['preInit'].pop()();
  }
}
// shouldRunNow refers to calling main(), not run().
var shouldRunNow = true;
if (Module['noInitialRun']) {
  shouldRunNow = false;
}
run();
// {{POST_RUN_ADDITIONS}}
  // {{MODULE_ADDITIONS}}