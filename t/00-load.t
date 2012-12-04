#!perl -T
use 5.006;
use strict;
use warnings FATAL => 'all';
use Test::More;

plan tests => 1;

BEGIN {
    use_ok( 'misc_perl_NY' ) || print "Bail out!\n";
}

diag( "Testing misc_perl_NY $misc_perl_NY::VERSION, Perl $], $^X" );
