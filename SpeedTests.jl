const max_iters = 20000000

fGlobal(i::Uint64) = i^2 + 2i + 1

type FunContainer
  f::Function
end

function testGlobal()
 for i=0:max_iters
   if fGlobal(uint64(i)) == 0
     error("Got 0")
   end
 end
end
# @elapsed testGlobal()
# 0.15, 0.13, 0.13

function testGlobalInContainer()
 const fc = FunContainer(fGlobal)
 for i=0:max_iters
   if fc.f(uint64(i)) == 0
     error("Got 0")
   end
 end
end
# @elapsed testGlobalInContainer()
# 2.01, 1.98

function testGlobalExtractedLocallyFromContainer()
 fc = FunContainer(fGlobal)
 const f = fc.f
 for i=0:max_iters
   if f(uint64(i)) == 0
     error("Got 0")
   end
 end
end
# @elapsed testGlobalExtractedLocallyFromContainer()
# 1.93, 1.92

function testLocalLambda()
 const f = i::Uint64 -> i^2 + 2i + 1
 fc = FunContainer(f)
 for i=0:max_iters
   if fc.f(uint64(i)) == 0
     error("Got 0")
   end
 end
end
# @elapsed testLocalLambda()
# 3.11, 3.07, 3.07

function testLocalFun()
 f(i::Uint64) = i^2 + 2i + 1
 const fc = FunContainer(f)
 for i=0:max_iters
   if fc.f(uint64(i)) == 0
     error("Got 0")
   end
 end
end
# @elapsed testLocalFun()
# 2.02, 1.97, 2.0

function testLocalFunNoContainer()
 f(i::Uint64) = i^2 + 2i + 1
 for i=0:max_iters
   if f(uint64(i)) == 0
     error("Got 0")
   end
 end
end
# @elapsed testLocalFunNoContainer()
# 1.91,  1.90, 1.92


