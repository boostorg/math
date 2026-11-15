// Copyright 2026 Nick Thompson
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// The exports that throw, which the bindings name with throws() and entry(), are wrapped: a C++ exception of this
// module becomes a RangeError with its message, anything else passes through unchanged, and either way the stack
// pointer goes back to where it was when this call came in. A WebAssembly exception leaves the stack where it was
// thrown, and a call made from a callback has to put it back to where its own caller left it.
{
  const fail = (sp, e) => {
    stackRestore(sp);
    if (!(e instanceof WebAssembly.Exception && e.is(getCppExceptionTag()))) throw e;
    let message;
    try { message = getExceptionMessage(e)[1]; } finally { decrementExceptionRefcount(e); }
    throw new RangeError(message);
  };
  const wrap = (f) => Object.assign(function () {
    const sp = stackSave();
    try { return f.apply(this, arguments); } catch (e) { fail(sp, e); }
  }, f);
  for (const name of new Set(Module.throwing)) {
    const [first, second] = name.split(/[ .]/);
    if (second === undefined) Module[first] = wrap(Module[first]);                       // f
    else if (first !== 'new') Module[first].prototype[second] = wrap(Module[first].prototype[second]);   // C.method
    else Module[second] = new Proxy(Module[second], {                                    // new C
      construct(target, args, newTarget) {
        const sp = stackSave();
        try { return Reflect.construct(target, args, newTarget); } catch (e) { fail(sp, e); }
      },
    });
  }
  delete Module.throwing;
}
