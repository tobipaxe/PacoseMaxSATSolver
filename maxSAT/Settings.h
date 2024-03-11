/********************************************************************************************
Copyright (c) 2017-2020, Tobias Paxian

dPermission is hereby granted, free of charge, to any person obtaining a copy of
    this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
********************************************************************************************/
 
#ifndef SETTINGS_H
#define SETTINGS_H

#include <stdint.h>
#include <iostream>
#include <string>

namespace Pacose {

enum EncodingType
{
  WARNERS,            // adder warners 1996, [Warners 1998?]
  BAILLEUX,           // totalizer [Bailleux & Boufkhad 2003]
  ASIN,               // [Asin et. al 2011] Robert Asin, Robert Nieuwenhuis, Albert Oliveras,
                      // Enric Rodriguez-Carbonell "Cardinality Networks: a theoretical and
                      // empirical study"
  OGAWA,              // modulo Totalizer [Ogawa et. al 2013]
  BAILLEUXW2,         // BailleuxW2
  WMTO,               // Weighted MaxSAT Totalizer
  MRWTO,              // Mixed Radix Weighted Totalizer
  MRWTO2,             // Mixed Radix Weighted Totalizer
  MRWTO19,            // Mixed Radix Weighted Totalizer 19 competition version
  MRWTO19_2,          // Mixed Radix Weighted Totalizer 19 competition version
  DGPW18,             // Dynamic Global Polynomial Watchdog [Paxian & Reimer 2018]
  HEURISTICQMAXSAT17, // auto heuristic [Koshi, 2014]
                      // selects between "warn, bail, ogaw"
  HEURISTIC20,        // HEURISTIC 20 Competition, choosing between warners and dgpw
  HEURISTIC1819,      // HEURISTIC 18 Competition, choosing between warners and
                      // dgpw
  QMAXSAT19
};

enum SolverType
{
  SOLVER,
  SIMPSOLVER,
  PARALLELSOLVER,
  MAXSATSOLVER
};

enum FormulaType
{
  NORMALSAT,
  MAXSAT,
  WEIGHTEDMAXSAT
};

enum EncodeStrategy
{
  ENCODEALL,
  ENCODEONLYIFNEEDED
};

enum PartitionStrategy
{
  NOPARTITION,
  GROUPBYWEIGHTADDATLAST,
  GROUPBYWEIGHT,
  GROUPBYBIGGESTREPEATINGENTRY,
  GROUPBYCOLUMNS
};

enum MultipleCascadeSolvingState
{
  SINGLECASCADE,
  MAINCASCADE,
  TARECASCADE,
  TARETARES
};

enum MultipleCascadeDivideStrategy
{
  SOLVEINNORMALCASCADEMODE,
  SORTEDNUMBEROFSOFTCLAUSES,
  RANDOMNUMBEROFSOFTCLAUSES,
  SOFTCLAUSESINORDER,
  SORTEDGOODDIVIDEPOINTS
};

enum InterimResult
{
  NOINTERIMRESULT,
  CUTATTOP,
  CUTATBOTTOM,
  CUTBOTH,
};

enum StructureInfo
{
  ISSAT,
  CONVERTTOMAXSAT,
  ISMAXSAT,
  DIVWEIGHTSBYDIVISOR,
  ISWEIGHTEDMAXSAT,
};

enum SorterType
{
  BITONIC,
  ODDEVEN,
  TOTALIZER,
};

enum DivideDGPWStrategy
{
  NODIVISION,
  USEONLYGCD,
  DIVIDEALL,
  DIVIDEALLSOLVEONLYWATCHDOGS,
  DIVIDEALLSOLVEGREEDY,
};

struct Settings
{
public:
  Settings()
       : //_encoding(HEURISTIC20),
        _encoding(DGPW18),
        _solverType(SOLVER),
        _formulaType(WEIGHTEDMAXSAT),
        _compression(-1),
        _printModel(true),
        _analyzeFormula(true),
        verbosity(0),
        cpuLimit(0),
        memLimit(0),
        maxWidth(0),
        networkType(TOTALIZER),
        preciseTarget(false),
        targetOpt(-1),
        analyze(true),
        encode01(false),
        lastPos1(true),
        base(2),
        groupHeuristic(1),
        percentOff(100),
        percentOffReinsert(true),
        equalWeight(0),
        partitionStrategy(GROUPBYWEIGHT),
        solveAtFirst(true),
        encodeStrategy(ENCODEONLYIFNEEDED),
        createGraphFile(""),
        onlyByTares(false),
        solveAsWeighted(true),
        formulaIsDivided(false),
        adderCaching(true),
        coneOfInfluence(true),
        exactBounding(true),
        plainVariant(false),
        dGPW(true),

        // multiple cascade mode
        mcDivideStrategy(SOLVEINNORMALCASCADEMODE),
        cascadeDivider(0),
        maxBucketSize(0),
        nOfCasc(0),
        tareCascadeOnlyByTares(false),
        sepHiWeight(false),
        weightPlusOne(false),
        exhaustiveWeightDistCheck(false),

        // interim results
        interimResult(NOINTERIMRESULT),
        divideDGPW(USEONLYGCD),

        minSize(0),
        partitionFactor(1),
        maxCnfFile(""),
        createSimplifiedWCNF(false),

        incremental(false),
        reuseDGPW(false),
        // greedyPrepro(1),
        greedyPrepro(0),
        greedyPPTimeLimit(600),
        greedyPPPropagationPerSecond(500000),
        greedyPPTimeoutFactor(18),
        greedyPPSSCSwitchFactor(75),
        greedyPPFixSCs(-1),
        greedyPPSATPercentage(65),
        greedyPPUNSATPercentage(65),
        greedyMinSizeOfSet(1),

        useGreedyPreInBetween(false),
        atLeastnEqualWeights(1),
        divisionMode(0),
        testIfDividable(2),
        divCheck(false),
        checkIfSolutionIsUnique(false),
        calculateAllSoftClauseCombinations(false),
        calculateAllSolutions(false),
        // CMS configuration Option
        reconf(99),
        simplify(false),
        useMaxPre2(false),
        maxPre2TimeOut(300),
        maxPre2Techniques("[bu]#[buvsrgcHVRTG]"),
        featureTest(0),
        proofFile("pacose_proof.pbp")
  {
    ResetDGPW();
    // ResetCore();
    // std::cout << "CARD=" << card << std::endl;
    verbosity = 0;
  }

  ~Settings() {}

  /**
   * @brief ResetCore
   *          sets all settings to standard values.
   */
  void ResetCore()
  {
    _encoding = HEURISTIC20;
    _solverType = SOLVER;
    _formulaType = WEIGHTEDMAXSAT;
    _compression = -1;
    _printModel = true;
    _analyzeFormula = true;
  }

  std::string ReturnEncodingString(EncodingType encoding)
  {
    if (encoding == WARNERS)
      return "warners";
    if (encoding == BAILLEUX)
      return "totalizer";
    if (encoding == ASIN)
      return "asin";
    if (encoding == OGAWA)
      return "ogawa";
    if (encoding == BAILLEUXW2)
      return "BailleuxW2";
    if (encoding == WMTO)
      return "WMTO";
    if (encoding == MRWTO)
      return "MRWTO";
    if (encoding == MRWTO2)
      return "MRWTO2";
    if (encoding == MRWTO19)
      return "MRWTO19";
    if (encoding == MRWTO19_2)
      return "MRWTO19_2";
    if (encoding == DGPW18)
      return "DGPW18";
    if (encoding == HEURISTICQMAXSAT17)
      return "H-QMSAT17";
    if (encoding == HEURISTIC20)
      return "H-Pac20";
    if (encoding == HEURISTIC1819)
      return "H-Pac1819";
    if (encoding == QMAXSAT19)
      return "QMSat19";
    return "invalid encoding";
  }

  struct CurrentCascade
  {
    bool _onlyWithAssumptions = false;
    bool _solveTares = true;
    uint32_t iteration = 0;
  } currentCascade;

  void ResetDGPW(void)
  {
    reconf = 99;
    // MaxSAT Settings
    _encoding = HEURISTIC20;
    _compression = -1;
    divisionMode = 0;
    greedyPPFixSCs = -1;
    simplify = false;
    networkType = TOTALIZER;
    preciseTarget = false;
    targetOpt = -1;
    greedyPrepro = 1;
    divCheck = false;
    greedyMinSizeOfSet = 1;
    useGreedyPreInBetween = false;
    atLeastnEqualWeights = 1;
    greedyPPTimeLimit = 600;
    greedyPPSATPercentage = 65;
    greedyPPUNSATPercentage = 65;
    greedyPPPropagationPerSecond = 500000;
    greedyPPTimeoutFactor = 18;
    greedyPPSSCSwitchFactor = 75;
    createSimplifiedWCNF = false;

    // incremental MaxSAT
    incremental = false;
    reuseDGPW = false;

    // Weighted MaxSAT Settings
    encode01 = false;
    lastPos1 = true;
    base = 2;
    groupHeuristic = 1;
    percentOff = 100;
    percentOffReinsert = true;
    equalWeight = 0;
    partitionStrategy = GROUPBYWEIGHT;
    solveAtFirst = true;
    encodeStrategy = ENCODEONLYIFNEEDED;
    createGraphFile = "";
    onlyByTares = false;
    solveAsWeighted = true;
    formulaIsDivided = false;
    partitionFactor = 1;
    testIfDividable = 0;
    checkIfSolutionIsUnique = false;
    calculateAllSoftClauseCombinations = false;
    calculateAllSolutions = false;

    // multiple cascade mode
    mcDivideStrategy = SOLVEINNORMALCASCADEMODE;
    cascadeDivider = 0;
    maxBucketSize = 0;
    nOfCasc = 0;
    tareCascadeOnlyByTares = false;
    sepHiWeight = false;
    weightPlusOne = true;
    minSize = 0;

    // interim results
    interimResult = NOINTERIMRESULT;

    featureTest = 0;
    analyze = true;

    // Paper Options
    adderCaching = true;
    coneOfInfluence = true;
    exactBounding = true;
    plainVariant = false;
    dGPW = true;

    SetPaperOptions();
  }

  void SetPaperOptions()
  {
    if (plainVariant || adderCaching || coneOfInfluence || exactBounding)
    {
      encodeStrategy = ENCODEALL;
      // partitionStrategy = NOPARTITION;
      lastPos1 = false;
      weightPlusOne = false;
      encode01 = false;
      solveAtFirst = true;
      solveAsWeighted = true;
      analyze = true;
      if (adderCaching && coneOfInfluence && exactBounding)
      {
        dGPW = true;
        plainVariant = false;
      }
    }
    if (adderCaching)
    {
      // partitionStrategy = GROUPBYWEIGHT;
      partitionStrategy = GROUPBYBIGGESTREPEATINGENTRY;
      // partitionStrategy = NOPARTITION;
      groupHeuristic = 1;
                    // partitionStrategy = GROUPBYCOLUMNS;
      //              createGraphFile = "graph";
    }
    if (coneOfInfluence)
    {
      encodeStrategy = ENCODEONLYIFNEEDED;
      lastPos1 = true;
    }
    if (exactBounding)
    {
      weightPlusOne = true;
      // had wrong value in competition version
      // slightly different from what we explained in the paper
      //              onlyByTares = true;
    }
  }

  void SetSolverType(SolverType solverType) { _solverType = solverType; }

  void SetEncoding(EncodingType encoding) { _encoding = encoding; }

  void SetFormulaType(FormulaType formulaType) { _formulaType = formulaType; }

  void SetCompression(int compression) { _compression = compression; }

  void SetPrintModel(bool printModel) { _printModel = printModel; }

  bool SetAnalyzeFormula() { return _analyzeFormula; }

  SolverType GetSolverType() { return _solverType; }

  EncodingType GetEncoding() { return _encoding; }

  FormulaType GetFormulaType() { return _formulaType; }

  int GetCompression() { return _compression; }

  bool GetPrintModel() { return _printModel; }

  bool GetAnalyzeFormula() { return _analyzeFormula; }

  void Print() const
  {
    if (verbosity == 0)
      return;

    std::cout << "c encoding...............: ";
    if (_encoding == HEURISTIC20)
    {
      std::cout << "HEURISTIC20" << std::endl;
    }
    else if (_encoding == DGPW18)
    {
      std::cout << "DGPW18" << std::endl;
    }
    else
    {
      std::cout << "other" << std::endl;
    }

    std::cout << "c solver type............: ";
    if (_solverType == SOLVER)
    {
      std::cout << "SOLVER" << std::endl;
    }
    else
    {
      std::cout << "other" << std::endl;
    }

    std::cout << "c formulaType............: ";
    if (_formulaType == WEIGHTEDMAXSAT)
    {
      std::cout << "WMS" << std::endl;
    }
    else
    {
      std::cout << "other" << std::endl;
    }

    std::cout << "c compression............: " << _compression << std::endl;

    std::cout << "c printModel.............: " << _printModel << std::endl;

    std::cout << "c analyzeFormula.........: " << _analyzeFormula << std::endl;

    std::cout << "c verbosity..............: " << verbosity << std::endl;
    std::cout << "c cpuLimit...............: " << cpuLimit << std::endl;
    std::cout << "c memLimit...............: " << memLimit << std::endl;
    std::cout << "c maxWidth...............: " << maxWidth << std::endl;
    std::cout << "c networkType............: ";
    if (networkType == TOTALIZER)
    {
      std::cout << "TOTALIZER" << std::endl;
    }
    else
    {
      std::cout << "other" << std::endl;
    }

    std::cout << "c preciseTarget..........: " << preciseTarget << std::endl;
    std::cout << "c targetOpt..............: " << targetOpt << std::endl;
    std::cout << "c analyze................: " << analyze << std::endl;
    std::cout << "c encode01...............: " << encode01 << std::endl;
    std::cout << "c lastPos1...............: " << lastPos1 << std::endl;
    std::cout << "c base...................: " << base << std::endl;
    std::cout << "c percentOff.............: " << percentOff << std::endl;
    std::cout << "c percentOffReinsert.....: " << percentOffReinsert << std::endl;
    std::cout << "c equalWeight............: " << equalWeight << std::endl;
    std::cout << "c partitionStrategy......: ";
    if (partitionStrategy == NOPARTITION)
    {
      std::cout << "NOPARTITION" << std::endl;
    }
    else if (partitionStrategy == GROUPBYWEIGHT)
    {
      std::cout << "GROUPBYWEIGHT" << std::endl;
    }
    else if (partitionStrategy == GROUPBYCOLUMNS)
    {
      std::cout << "GROUPBYCOLUMNS" << std::endl;
    }
        else if (partitionStrategy == GROUPBYBIGGESTREPEATINGENTRY)
    {
      std::cout << "GROUPBYBIGGESTREPEATINGENTRY" << std::endl;
    }
    else
    {
      std::cout << "other" << std::endl;
    }
    std::cout << "c groupHeuristic.........: " << groupHeuristic << std::endl;
    std::cout << "c solveAtFirst...........: " << solveAtFirst << std::endl;
    std::cout << "c encodeStrategy.........: ";
    if (encodeStrategy == ENCODEONLYIFNEEDED)
    {
      std::cout << "ENCODEONLYIFNEEDED" << std::endl;
    }
    else
    {
      std::cout << "other" << std::endl;
    }
    std::cout << "c createGraphFile........: " << createGraphFile << std::endl;
    std::cout << "c onlyByTares............: " << onlyByTares << std::endl;
    std::cout << "c solveAsWeighted........: " << solveAsWeighted << std::endl;
    std::cout << "c formulaIsDivided.......: " << formulaIsDivided << std::endl;
    std::cout << "c adderCaching...........: " << adderCaching << std::endl;
    std::cout << "c coneOfInfluence........: " << coneOfInfluence << std::endl;
    std::cout << "c exactBounding..........: " << exactBounding << std::endl;
    std::cout << "c plainVariant...........: " << plainVariant << std::endl;
    std::cout << "c dGPW...................: " << dGPW << std::endl;

    // multiple cascade mode
    std::cout << "c mcDivideStrategy.......: ";
    if (mcDivideStrategy == SOLVEINNORMALCASCADEMODE)
    {
      std::cout << "SOLVEINNORMALCASCADEMODE" << std::endl;
    }
    else
    {
      std::cout << "other" << std::endl;
    }
    std::cout << "c cascadeDivider.........: " << cascadeDivider << std::endl;
    std::cout << "c maxBucketSize..........: " << maxBucketSize << std::endl;
    std::cout << "c nOfCasc................: " << nOfCasc << std::endl;
    std::cout << "c tareCascadeOnlyByTares.: " << tareCascadeOnlyByTares << std::endl;
    std::cout << "c sepHiWeight............: " << sepHiWeight << std::endl;
    std::cout << "c weightPlusOne..........: " << weightPlusOne << std::endl;

    // interim results
    std::cout << "c interimResult..........: ";
    if (interimResult == NOINTERIMRESULT)
    {
      std::cout << "NOINTERIMRESULT" << std::endl;
    }
    else
    {
      std::cout << "other" << std::endl;
    }
    
    std::cout << "c divideDGPW.............: ";
    if (divideDGPW == NODIVISION)
    {
      std::cout << "NODIVISION" << std::endl;
    }
    else
    {
      std::cout << "other" << std::endl;
    }
    
    std::cout << "c minSize................: " << minSize << std::endl;
    std::cout << "c partitionFactor........: " << partitionFactor << std::endl;
    std::cout << "c maxCnfFile.............: " << maxCnfFile << std::endl;
    std::cout << "c createSimplifiedWCNF...: " << createSimplifiedWCNF << std::endl;
    std::cout << "c incremental............: " << incremental << std::endl;

    std::cout << "c reuseDGPW..............: " << reuseDGPW << std::endl;
    std::cout << "c greedyPrepro...........: " << greedyPrepro << std::endl;
    std::cout << "c greedyPPTimeLimit......: " << greedyPPTimeLimit << std::endl;
    std::cout << "c greedyPPPropagationPerS: " << greedyPPPropagationPerSecond << std::endl;
    std::cout << "c greedyPPTimeoutFactor..: " << greedyPPTimeoutFactor << std::endl;
    std::cout << "c greedyPPSSCSwitchFactor: " << greedyPPSSCSwitchFactor << std::endl;
    std::cout << "c greedyPPFixSCs.........: " << greedyPPFixSCs << std::endl;
    std::cout << "c greedyPPSATPercentage..: " << greedyPPSATPercentage << std::endl;
    std::cout << "c greedyPPUNSATPercentage: " << greedyPPUNSATPercentage << std::endl;
    std::cout << "c greedyMinSizeOfSet.....: " << greedyMinSizeOfSet << std::endl;
    std::cout << "c useGreedyPreInBetween..: " << useGreedyPreInBetween << std::endl;
    std::cout << "c atLeastnEqualWeights...: " << atLeastnEqualWeights << std::endl;
    std::cout << "c divisionMode...........: " << divisionMode << std::endl;
    std::cout << "c testIfDividable........: " << testIfDividable << std::endl;
    std::cout << "c divCheck...............: " << divCheck << std::endl;
    std::cout << "c checkIfSolutionIsUnique: " << checkIfSolutionIsUnique << std::endl;
    std::cout << "c calcAll SC combinations: " << calculateAllSoftClauseCombinations << std::endl;
    std::cout << "c calculateAllSolutions..: " << calculateAllSolutions << std::endl;
    std::cout << "c reconf.................: " << reconf << std::endl;
    // CMS configuration Option
    std::cout << "c simplify...............: " << simplify << std::endl;
    std::cout << "c featureTest............: " << featureTest << std::endl;

    // std::cout << "c Plain Variant..........: "
    //         << (plainVariant ? "true" : "false") << std::endl;
    // std::cout << "c Adder Caching..........: "
    //           << (adderCaching ? "true" : "false") << std::endl;
    // std::cout << "c Cone Of Influence......: "
    //           << (coneOfInfluence ? "true" : "false") << std::endl;
    // std::cout << "c Exact Bounding.........: "
    //           << (exactBounding ? "true" : "false") << std::endl;
    // std::cout << "c DGPW combined..........: " << (dGPW ? "true" : "false")
    //           << std::endl;
    // std::cout << "c verbosity..............: " << verbosity << std::endl;
    // std::cout << "c ------------------------" << std::endl;
    // std::cout << "c MaxSAT options:" << std::endl;
    // std::cout << "c weight plus 1..........: " << weightPlusOne << std::endl;
    // std::cout << "c network type...........: ";
    // switch (networkType)
    // {
    // case BITONIC:
    //   std::cout << "Bitonic sorter" << std::endl;
    //   break;
    // case ODDEVEN:
    //   std::cout << "Odd-Even sorter" << std::endl;
    //   break;
    // case TOTALIZER:
    //   std::cout << "Totalizer" << std::endl;
    //   break;
    // }
    // //          std::cout << "c decision strategies....: ";
    // std::cout << "c encode 01..............: " << (encode01 ? "true" : "false")
    //           << std::endl;
    // //          std::cout << "c sort soft clauses......: ";
    // std::cout << "c base...................: " << base << std::endl;
    // std::cout << "c analyze................: " << (analyze ? "true" : "false")
    //           << std::endl;
    // std::cout << "c partitionStrategy......: " << partitionStrategy
    //           << std::endl;
    // std::cout << "c   heuristic............: " << groupHeuristic << std::endl;

    // std::cout << "c defined target optimum.: " << targetOpt << std::endl;
  }

  //  private:
  EncodingType _encoding;
  SolverType _solverType;
  FormulaType _formulaType;
  // compression of encoding due to comp in QMaxSAT;
  // 0 highest compression
  int _compression;
  bool _printModel;
  bool _analyzeFormula;

  // General Settings
  uint32_t verbosity;
  double cpuLimit;
  uint32_t memLimit;
  uint32_t maxWidth;
  SorterType networkType;
  bool preciseTarget;
  int32_t targetOpt;

  // Weighted MaxSAT Settings
  bool analyze;
  bool encode01;
  bool lastPos1;
  uint32_t base;
  uint32_t groupHeuristic;
  uint32_t percentOff;
  bool percentOffReinsert;
  uint32_t equalWeight;
  PartitionStrategy partitionStrategy;
  bool solveAtFirst;
  EncodeStrategy encodeStrategy;
  std::string createGraphFile;
  bool onlyByTares;
  bool solveAsWeighted;
  bool formulaIsDivided;

  // Paper Options
  bool adderCaching;
  bool coneOfInfluence;
  bool exactBounding;
  bool plainVariant;
  bool dGPW;

  // multiple cascade mode
  MultipleCascadeDivideStrategy mcDivideStrategy;
  uint32_t cascadeDivider;
  uint32_t maxBucketSize;
  uint32_t nOfCasc;
  bool tareCascadeOnlyByTares;
  bool sepHiWeight;
  bool weightPlusOne;
  bool exhaustiveWeightDistCheck;

  // interim results
  InterimResult interimResult;
  DivideDGPWStrategy divideDGPW;

  uint32_t minSize;
  uint32_t partitionFactor;
  std::string maxCnfFile;
  //  bool useGreedyPrepro;
  bool createSimplifiedWCNF;

  bool incremental;
  bool reuseDGPW;
  int greedyPrepro;
  int greedyPPTimeLimit;
  int greedyPPPropagationPerSecond;
  int greedyPPTimeoutFactor;
  int greedyPPSSCSwitchFactor;
  int greedyPPFixSCs;
  int greedyPPSATPercentage;
  int greedyPPUNSATPercentage;
  int greedyMinSizeOfSet;

  bool useGreedyPreInBetween;
  uint32_t atLeastnEqualWeights;
  int divisionMode;
  int testIfDividable;
  bool divCheck;
  bool checkIfSolutionIsUnique;
  bool calculateAllSoftClauseCombinations;
  bool calculateAllSolutions;

  // CMS configuration Option
  uint32_t reconf;
  bool simplify;

  // maxpre2
  bool useMaxPre2;
  double maxPre2TimeOut;
  std::string maxPre2Techniques;

  int featureTest;

  // proofing
  std::string proofFile;
};
} // Namespace Pacose

#endif // SETTINGS_H
