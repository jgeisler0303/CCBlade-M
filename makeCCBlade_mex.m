function makeCCBlade_mex
mex_name= 'CCBlade_mex';
mex_file= [mex_name '.cpp'];

v= ver;
is_matlab= ~strcmp(v(1).Name, 'Octave');
if is_matlab
    clear mex
else
    clear(mex_name);
end

options= {};
if is_matlab
    options= [options {'CXXFLAGS="$CXXFLAGS -std=c++11 -Wall -fdiagnostics-show-option"'}];
else
    old_cxxflags= getenv('CXXFLAGS');
    [status, cxxflags]= system('mkoctfile --print  CXXFLAGS');
    cxxflags= [cxxflags ' -std=c++11 -Wall -fdiagnostics-show-option '];
    cxxflags= strrep(cxxflags, char(10), ' ');
    setenv('CXXFLAGS', cxxflags);
end

mex(options{:}, mex_file);

if ~is_matlab
    setenv('CXXFLAGS', old_cxxflags);
end
