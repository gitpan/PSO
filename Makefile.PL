use 5.008008;
use ExtUtils::MakeMaker;
WriteMakefile(
    NAME              => 'PSO',
    VERSION_FROM      => 'lib/PSO.pm',
    PREREQ_PM         => {},
    ($] >= 5.005 ?
      (ABSTRACT_FROM  => 'lib/PSO.pm',
       AUTHOR         => 'W. Kyle Schlansker <kylesch@gmail.com>') : ()),
);
