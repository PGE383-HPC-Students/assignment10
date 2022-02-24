#!/usr/bin/env julia

# Copyright 2022 John T. Foster
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import assignment10
using Test

@testset "assignment10.jl" begin
    params = [2500.0, 1.0e6, 14.504, 518.67, 2000, 0.95, 2.0, 679.66, 0.056, 1.e2]
    sol = assignment10.gas_solver(params...)
    @test all(isapprox.(getindex.(sol.u, 1)[10:15],[2312.059, 2267.440, 2229.440, 2169.752, 2129.861, 2096.855], atol=1.0e-3))
end
