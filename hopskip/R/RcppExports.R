# This file was generated by Rcpp::compileAttributes
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

intersections <- function(sunitsx, sunitsy, sendpointsx, sendpointsy, sstreetsp0, sstreetsp1) {
    .Call('hopskip_intersections', PACKAGE = 'hopskip', sunitsx, sunitsy, sendpointsx, sendpointsy, sstreetsp0, sstreetsp1)
}

simple_hazard <- function(pairwise_distanceS, street_matrixS, parametersS, callback) {
    .Call('hopskip_simple_hazard', PACKAGE = 'hopskip', pairwise_distanceS, street_matrixS, parametersS, callback)
}

bugs <- function(pairwise_distanceS, parametersS) {
    .Call('hopskip_bugs', PACKAGE = 'hopskip', pairwise_distanceS, parametersS)
}

TestExcessGrowthDistribution <- function() {
    .Call('hopskip_TestExcessGrowthDistribution', PACKAGE = 'hopskip')
}

TestCallback <- function(callback) {
    .Call('hopskip_TestCallback', PACKAGE = 'hopskip', callback)
}

hilbertXY2D <- function(sx, sy, sn) {
    .Call('hopskip_hilbertXY2D', PACKAGE = 'hopskip', sx, sy, sn)
}

