use Test::More tests => 4;
BEGIN { use_ok('PSO') };

my %test_params = (
	numParticles   => 4,
	numNeighbors   => 3,
	maxIterations  => 1000,
	dimensions     => 4,
	deltaMin       => -2.0,
	deltaMax       =>  4.0,
	meWeight       => 2.0,
	meMin          => 0.0,
	meMax          => 1.0,
	themWeight     => 2.0,
	themMin        => 0.0,
	themMax        => 1.0,
	exitFitness    => 1.0,
	verbose        => 1,
);

# simple test function to sum the position values up to 3.5
my $testValue = 3.5;
sub test_fitness_function(@) {
        my (@arr) = (@_);
        my $sum = 0;
	my $ret = 0;
        foreach my $val (@arr) {
                $sum += $val;
        }
	if($sum <= $testValue + 0.01 && $sum >= $testValue - 0.01) {
		$ret = 1;
	}
	else {
		$ret = 0;
	}
	return $ret;
}


ok( pso_set_params(\%test_params) == 0 );
ok( pso_register_fitness_function('test_fitness_function') == 0);
ok( pso_optimize() == 0 );
