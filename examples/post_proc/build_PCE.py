from operator import itemgetter

import numpy as np
import openturns as ot


def construct_PCE_ot(
    training_input, training_output, marginals, copula, degree, LARS=True
):
    ##########################INPUTS##########################
    ##########################################################
    Nt = len(training_input)
    if len(training_input.shape) > 1:
        Nvar = training_input.shape[1]
    else:
        Nvar = 1

    # Define Sample
    outputSample = ot.Sample(Nt, 1)
    for i in range(Nt):
        outputSample[i, 0] = training_output[i]

    # Define Collection and PDFs
    polyColl = ot.PolynomialFamilyCollection(Nvar)
    collection = ot.DistributionCollection(Nvar)
    marginal = {}
    UncorrelatedInputSample = ot.Sample(Nt, Nvar)

    if Nvar > 1:
        for i in range(Nvar):
            varSample = ot.Sample(Nt, 1)
            for j in range(Nt):
                varSample[j, 0] = training_input[j, i]
                UncorrelatedInputSample[j, i] = training_input[j, i]
            minValue = varSample.getMin()[0]
            maxValue = varSample.getMax()[0]
            if marginals[i] == "gaussian" or marginals[i] == "normal":
                marginal[i] = ot.NormalFactory().build(varSample)
            elif marginals[i] == "uniform":
                marginal[i] = ot.Uniform(
                    minValue - minValue / 100.0, maxValue + maxValue / 100.0
                )
            elif marginals[i] == "kernel":
                marginal[i] = ot.KernelSmoothing().build(varSample)
            else:
                print(
                    "WARNING: couldn't find distribution '"
                    + str(marginals[i])
                    + "', applied kernel smoothing instead"
                )
                marginal[i] = ot.KernelSmoothing().build(varSample)

            collection[i] = ot.Distribution(marginal[i])
    else:
        varSample = ot.Sample(Nt, 1)
        for j in range(Nt):
            varSample[j, 0] = training_input[j]
            UncorrelatedInputSample[j, 0] = training_input[j]
        minValue = varSample.getMin()[0]
        maxValue = varSample.getMax()[0]
        if marginals[i] == "gaussian" or marginals[i] == "normal":
            marginal[0] = ot.NormalFactory().build(varSample)
        elif marginals[i] == "uniform":
            marginal[0] = ot.Uniform(
                minValue - minValue / 100.0, maxValue + maxValue / 100.0
            )
        elif marginals[i] == "kernel":
            marginal[0] = ot.KernelSmoothing().build(varSample)
        else:
            print(
                "WARNING: couldn't find distribution '"
                + str(marginals[i])
                + "', applied kernel smoothing instead"
            )
            marginal[0] = ot.KernelSmoothing().build(varSample)
        collection[0] = ot.Distribution(marginal[0])

    if copula == "independent":
        copula = ot.IndependentCopula(Nvar)
    elif copula == "gaussian" or copula == "normal":
        inputSample = ot.Sample(training_input)
        copula = ot.NormalCopulaFactory().build(inputSample)
    else:
        print(
            "WARNING: couldn't find copula '"
            + str(copula)
            + "', applied independent copula instead"
        )
        copula = ot.IndependentCopula(Nvar)

    # UncorrelatedInputDistribution = ot.ComposedDistribution(collection,ot.Copula(copula))
    UncorrelatedInputDistribution = ot.ComposedDistribution(collection, copula)

    # Calcul des polynomes du chaos
    for v in range(0, Nvar):
        marginalv = UncorrelatedInputDistribution.getMarginal(v)
        if marginals[i] == "kernel":
            # Works with arbitrary PDF
            basisAlgorithm = ot.AdaptiveStieltjesAlgorithm(marginalv)
            polyColl[v] = ot.StandardDistributionPolynomialFactory(basisAlgorithm)
        else:
            # Works with standard PDF: gaussian, uniform, ..
            polyColl[v] = ot.StandardDistributionPolynomialFactory(marginalv)

    # Definition de la numerotation des coefficients des polynomes du chaos
    enumerateFunction = ot.LinearEnumerateFunction(Nvar)
    # enumerateFunction = HyperbolicAnisotropicEnumerateFunction(Nvar,0.4)
    # Creation de la base des polynomes multivaries en fonction de la numerotation
    #                     et des bases desiree
    multivariateBasis = ot.OrthogonalProductPolynomialFactory(
        polyColl, enumerateFunction
    )
    # Number of PC terms
    P = enumerateFunction.getStrataCumulatedCardinal(degree)
    # Troncature
    adaptativeStrategy = ot.FixedStrategy(multivariateBasis, P)

    if LARS:
        # Evaluation Strategy : LARS
        basisSequenceFactory = ot.LARS()
        fittingAlgorithm = ot.CorrectedLeaveOneOut()
        approximationAlgorithm = ot.LeastSquaresMetaModelSelectionFactory(
            basisSequenceFactory, fittingAlgorithm
        )

        # Approximation method for PCE coefficients
        projectionStrategy = ot.LeastSquaresStrategy(
            UncorrelatedInputSample, outputSample, approximationAlgorithm
        )
        #
        algo = ot.FunctionalChaosAlgorithm(
            UncorrelatedInputSample,
            outputSample,
            UncorrelatedInputDistribution,
            adaptativeStrategy,
            projectionStrategy,
        )
    else:
        #
        wei_exp = ot.MonteCarloExperiment(
            UncorrelatedInputDistribution, UncorrelatedInputSample.getSize()
        )
        X_UncorrelatedInputSample, weights = wei_exp.generateWithWeights()
        projectionStrategy = ot.LeastSquaresStrategy()
        algo = ot.FunctionalChaosAlgorithm(
            X_UncorrelatedInputSample,
            weights,
            outputSample,
            UncorrelatedInputDistribution,
            adaptativeStrategy,
            projectionStrategy,
        )

    algo.run()
    polynomialChaosResult = algo.getResult()
    metamodel = polynomialChaosResult.getMetaModel()
    enumerateFunction = enumerateFunction

    return polynomialChaosResult
