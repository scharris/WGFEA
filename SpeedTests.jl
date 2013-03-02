fGlobal(i::Uint64) = i^2 + 2i + 1

type FunContainer
  f::Function
end

function testGlobal()
 sum::Uint64 = 0
 for i=0:20000000
   if fGlobal(uint64(i)) == 0
     error("Got 0")
   else
    sum += i
   end
 end
 sum
end
# @elapsed testGlobal()
# 0.15, 0.13, 0.13

function testGlobalInContainer()
 sum::Uint64 = 0
 const fc = FunContainer(fGlobal)
 for i=0:20000000
   if fc.f(uint64(i)) == 0
     error("Got 0")
   else
    sum += i
   end
 end
 sum
end
# @elapsed testGlobalInContainer()
# 2.01, 1.98

function testLocalLambda()
 sum::Uint64 = 0
 const f = i::Uint64 -> i^2 + 2i + 1
 fc = FunContainer(f)
 for i=0:20000000
   if fc.f(uint64(i)) == 0
     error("Got 0")
   else
    sum += i
   end
 end
 sum
end
# @elapsed testLocalLambda()
# 3.11, 3.07, 3.07

function testLocalFun()
 sum::Uint64 = 0
 f(i::Uint64) = i^2 + 2i + 1
 const fc = FunContainer(f)
 for i=0:20000000
   if fc.f(uint64(i)) == 0
     error("Got 0")
   else
    sum += i
   end
 end
 sum
end
# @elapsed testLocalFun()
# 2.02, 1.97, 2.0

