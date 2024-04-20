# Meander Factorization Toolkit

## Overview
This repository contains a collection of specialized functions and classes designed to facilitate the exploration of meanders and the reproduction of existing results within the field. It is important to note that this software is not a complete application. The library does not include a top-level application or API; it is up to the user to implement the integration into their projects.
The code herein is developed for the article "Prime Factorisation of Meanders" and implements the basic computational apparatus described in that work. The preprint of the paper is available on [arXiv](https://arxiv.org/abs/2112.10289).

## Prerequisites
The toolkit utilizes the Boost C++ Libraries, specifically the `boost/multiprecision` component, to handle large numbers required in calculations. Users will need to ensure that the Boost libraries are properly installed and configured in their development environment.

## Structure and Examples
The main examples of how to utilize this toolkit can be found in `source.cpp`. This file includes:
- An example of calculating irreducible meander numbers.
- An example of computing numbers for iterated snakes.

Both examples also present the output of these calculations, serving as a practical guide to what users can expect when leveraging this toolkit for similar tasks.