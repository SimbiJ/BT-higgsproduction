#!/usr/bin/env python3
from pySecDec.loop_integral import LoopIntegralFromPropagators, LoopPackage
from pySecDec.code_writer import sum_package, Coefficient


REPLACEMENT_RULES = [
    ('p*p', 0),
    ('q*q', 0),
    ('h*h', 's12'),
    ('p*q', 's12/2'),
    ('m', 'msq')
]



def generate_loops():
    li1 = LoopIntegralFromPropagators(
        loop_momenta = ['l'],
        external_momenta = ['p','q', 'h'],

        propagators = ['l**2-m**2', '(l-q)**2-m**2', '(l+p)**2-m**2', 'l**2'],
        powerlist = [1,1,1,-1],

        replacement_rules = REPLACEMENT_RULES,

        dimensionality = '4-2*eps'
    )


    li2 = LoopIntegralFromPropagators(
        loop_momenta = ['l'],
        external_momenta = ['p','q', 'h'],

        Lorentz_indices = ['mu', 'nu'],

        propagators = ['l**2-m**2', '(l-q)**2-m**2', '(l+p)**2-m**2'],
        powerlist = [1,1,1],

        numerator = '(l(mu)*p(mu))*(l(nu)*q(nu))',

        replacement_rules = REPLACEMENT_RULES,

        dimensionality = '4-2*eps'
    )


    li3 = LoopIntegralFromPropagators(
        loop_momenta = ['l'],
        external_momenta = ['p','q', 'h'],

        propagators = ['l**2-m**2', '(l-q)**2-m**2', '(l+p)**2-m**2'],
        powerlist = [1,1,1],

        replacement_rules = REPLACEMENT_RULES,

        dimensionality = '4-2*eps'
    )

    return li1, li2, li3



def sum_integrals(li1, li2, li3):
    sum_package(
        name = 'sigma_num',

        package_generators = [
            LoopPackage(name='num_l_sq', loop_integral=li1, requested_order=[0]),
            LoopPackage(name='num_l_lin', loop_integral=li2, requested_order=[0]),
            LoopPackage(name='num_const', loop_integral=li3, requested_order=[0])
        ],

        regulators = ['eps'],

        requested_orders = [0],

        real_parameters = ['s12', 'msq'],

        coefficients = [[
            Coefficient(['4*msq*(5-(4-2*eps))'], ['(4-2*eps)-2'], ['s12', 'msq']),
            Coefficient(['-32*msq'], ['s12*((4-2*eps)-2)'], ['s12', 'msq']),
            Coefficient(['4*msq*(((4-2*eps)-1)*msq**2+(2-(4-2*eps))*s12/2)'], ['(4-2*eps)-2'], ['s12', 'msq'])
        ]] # Coefficient(numerator, denominator, parameters)
    )



def main():
    li1, li2, li3 = generate_loops()
    sum_integrals(li1, li2, li3)



if __name__ == "__main__":
    main()
