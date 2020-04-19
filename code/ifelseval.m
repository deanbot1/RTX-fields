function V = ifelseval(test,trueval,falseval)
% little utility function that returns trueval is test is true, falseval otherwise
% it also turns unitary cell arrays to the thing that's in them.

if test
	V = trueval;
else
	V = falseval;
end

if iscell(V) & length(V)==1
	V = V{1};
end
