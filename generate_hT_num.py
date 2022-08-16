#!/usr/bin/env python3
from pySecDec.loop_integral import LoopIntegralFromPropagators, LoopPackage, loop_regions
from pySecDec.code_writer import sum_package, Coefficient



def define_integral():
    replacement_rules = [
        ('p*p', 0),
        ('q*q', 0),
        ('h*h', 's12'),
        ('p*q', 's12/2'),
        ('1/(m**2)', 'invmsq'),
    ]

    C_0 = LoopIntegralFromPropagators(

    loop_momenta = ['l'],
    external_momenta = ['p','q', 'h'],

    propagators = ['l**2-m**2', '(l-q)**2-m**2', '(l+p)**2-m**2'],
    powerlist = [1, 1, 1],

    replacement_rules = replacement_rules,

    regulators = ['eps'],

    dimensionality = '4-2*eps'

    )

    return C_0



def split_sum_regions(C_0):
    regions = loop_regions(

        name = 'heavyTop',

        loop_integral = C_0,

        smallness_parameter = 'invmsq',

        expansion_by_regions_order = 2,
    )


    sum_package(

        name = 'heavyTop',

        package_generators = regions,

        regulators = ['eps'],

        requested_orders = [0],

        real_parameters = ['s12', 'invmsq'],
    )



def main():
    C_0 = define_integral()
    split_sum_regions(C_0)



if __name__ == "__main__":
    main()
