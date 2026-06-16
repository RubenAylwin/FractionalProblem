#ifndef SOLVE_RL_PROBLEM_H
#define SOLVE_RL_PROBLEM_H

#include <string>
#include <MyTypes.h>

void solveRLProblem(int order,
                    const std::string &typeD,
                    const std::string &typeQ,
                    const BEM::CVector &dVec,
                    const BEM::CVector &qVec,
                    const BEM::CVector &rhsVec,
                    int ms);

#endif // SOLVE_RL_PROBLEM_H
