/*
 * Copyright (C) 2021 Peter G. Jensen <root@petergjoel.dk>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#define BOOST_TEST_MODULE UnorderedLoad

#include <boost/test/unit_test.hpp>

#include "SimpleTree.h"


const std::string simple_unordered_strategy = "{\"version\":1.0,\"type\":\"state->regressor\",\"representation\":\"map\",\"actions\":{"
"	\"0\":\"Controller._id4->Controller._id2 { 1, tau, minute_clock := TIME_OFFSET, initValues() }\","
"	\"1\":\"Controller._id2->Controller.Wait { 0 >= 10 && any_on(), tau, consumed_power := 0.625000 + 0.187500 * 0, heat_produced := consumed_power * calculateCOP(), setMassFlow() }\","
"	\"2\":\"Controller._id2->Controller.Wait { 1 >= 10 && any_on(), tau, consumed_power := 0.625000 + 0.187500 * 1, heat_produced := consumed_power * calculateCOP(), setMassFlow() }\","
"	\"3\":\"Controller._id2->Controller.Wait { 2 >= 10 && any_on(), tau, consumed_power := 0.625000 + 0.187500 * 2, heat_produced := consumed_power * calculateCOP(), setMassFlow() }\","
"	\"4\":\"Controller._id2->Controller.Wait { 3 >= 10 && any_on(), tau, consumed_power := 0.625000 + 0.187500 * 3, heat_produced := consumed_power * calculateCOP(), setMassFlow() }\","
"	\"5\":\"Controller._id2->Controller.Wait { 4 >= 10 && any_on(), tau, consumed_power := 0.625000 + 0.187500 * 4, heat_produced := consumed_power * calculateCOP(), setMassFlow() }\","
"	\"6\":\"Controller._id2->Controller.Wait { 5 >= 10 && any_on(), tau, consumed_power := 0.625000 + 0.187500 * 5, heat_produced := consumed_power * calculateCOP(), setMassFlow() }\","
"	\"7\":\"Controller._id2->Controller.Wait { 6 >= 10 && any_on(), tau, consumed_power := 0.625000 + 0.187500 * 6, heat_produced := consumed_power * calculateCOP(), setMassFlow() }\","
"	\"9\":\"Controller._id2->Controller.Wait { 8 >= 10 && any_on(), tau, consumed_power := 0.625000 + 0.187500 * 8, heat_produced := consumed_power * calculateCOP(), setMassFlow() }\","
"	\"10\":\"Controller._id2->Controller.Wait { 9 >= 10 && any_on(), tau, consumed_power := 0.625000 + 0.187500 * 9, heat_produced := consumed_power * calculateCOP(), setMassFlow() }\","
"	\"11\":\"Controller._id2->Controller.Wait { 10 >= 10 && any_on(), tau, consumed_power := 0.625000 + 0.187500 * 10, heat_produced := consumed_power * calculateCOP(), setMassFlow() }\","
"	\"12\":\"Controller._id2->Controller.Wait { 1, tau, consumed_power := 0, heat_produced := 0 }\","
"	\"13\":\"WAIT\""
"},\"statevars\":["
"	\"round(minute_clock) \""
"],\"pointvars\":["
"],\"locationnames\":{"
"	\"Fetch_Data.location\":{"
"		\"0\":\"_id5\""
"	},"
"	\"Room1.location\":{"
"		\"0\":\"_id0\""
"	},"
"	\"Room2.location\":{"
"		\"0\":\"_id0\""
"	},"
"	\"Room3.location\":{"
"		\"0\":\"_id0\""
"	},"
"	\"Room4.location\":{"
"		\"0\":\"_id0\""
"	},"
"	\"Optimization.location\":{"
"		\"0\":\"_id1\""
"	},"
"	\"Controller.location\":{"
"		\"0\":\"_id2\","
"		\"1\":\"Wait\","
"		\"2\":\"_id4\""
"	}"
"},\"regressors\":{"
"		\"(30)\":"
"			{\"type\":\"act->point->val\",\"representation\":\"simpletree\",\"minimize\":1,\"regressor\":"
"				{"
"				\"11\" : 1.294636921229531,"
"				\"12\" : 3.091107032150207"
"				}"
"			},"
"		\"(210)\":"
"			{\"type\":\"act->point->val\",\"representation\":\"simpletree\",\"minimize\":1,\"regressor\":"
"				{"
"				\"11\" : 1.009248500699941,"
"				\"12\" : 1.028772620412916"
"				}"
"			},"
"		\"(15)\":"
"			{\"type\":\"act->point->val\",\"representation\":\"simpletree\",\"minimize\":1,\"regressor\":"
"				{"
"				\"11\" : 2.855407148131672,"
"				\"12\" : 1.352930166302546"
"				}"
"			}"
"	}"
"}";




const std::string unordered_strategy = "{\"version\":1.0,\"type\":\"state->regressor\",\"representation\":\"map\",\"actions\":{"
"	\"0\":\"Controller._id4->Controller._id2 { 1, tau, minute_clock := TIME_OFFSET, initValues() }\","
"	\"1\":\"Controller._id2->Controller.Wait { 0 >= 10 && any_on(), tau, consumed_power := 0.625000 + 0.187500 * 0, heat_produced := consumed_power * calculateCOP(), setMassFlow() }\","
"	\"2\":\"Controller._id2->Controller.Wait { 1 >= 10 && any_on(), tau, consumed_power := 0.625000 + 0.187500 * 1, heat_produced := consumed_power * calculateCOP(), setMassFlow() }\","
"	\"3\":\"Controller._id2->Controller.Wait { 2 >= 10 && any_on(), tau, consumed_power := 0.625000 + 0.187500 * 2, heat_produced := consumed_power * calculateCOP(), setMassFlow() }\","
"	\"4\":\"Controller._id2->Controller.Wait { 3 >= 10 && any_on(), tau, consumed_power := 0.625000 + 0.187500 * 3, heat_produced := consumed_power * calculateCOP(), setMassFlow() }\","
"	\"5\":\"Controller._id2->Controller.Wait { 4 >= 10 && any_on(), tau, consumed_power := 0.625000 + 0.187500 * 4, heat_produced := consumed_power * calculateCOP(), setMassFlow() }\","
"	\"6\":\"Controller._id2->Controller.Wait { 5 >= 10 && any_on(), tau, consumed_power := 0.625000 + 0.187500 * 5, heat_produced := consumed_power * calculateCOP(), setMassFlow() }\","
"	\"7\":\"Controller._id2->Controller.Wait { 6 >= 10 && any_on(), tau, consumed_power := 0.625000 + 0.187500 * 6, heat_produced := consumed_power * calculateCOP(), setMassFlow() }\","
"	\"9\":\"Controller._id2->Controller.Wait { 8 >= 10 && any_on(), tau, consumed_power := 0.625000 + 0.187500 * 8, heat_produced := consumed_power * calculateCOP(), setMassFlow() }\","
"	\"10\":\"Controller._id2->Controller.Wait { 9 >= 10 && any_on(), tau, consumed_power := 0.625000 + 0.187500 * 9, heat_produced := consumed_power * calculateCOP(), setMassFlow() }\","
"	\"11\":\"Controller._id2->Controller.Wait { 10 >= 10 && any_on(), tau, consumed_power := 0.625000 + 0.187500 * 10, heat_produced := consumed_power * calculateCOP(), setMassFlow() }\","
"	\"12\":\"Controller._id2->Controller.Wait { 1, tau, consumed_power := 0, heat_produced := 0 }\","
"	\"13\":\"WAIT\""
"},\"statevars\":["
"	\"round(minute_clock) \""
"],\"pointvars\":["
"],\"locationnames\":{"
"	\"Fetch_Data.location\":{"
"		\"0\":\"_id5\""
"	},"
"	\"Room1.location\":{"
"		\"0\":\"_id0\""
"	},"
"	\"Room2.location\":{"
"		\"0\":\"_id0\""
"	},"
"	\"Room3.location\":{"
"		\"0\":\"_id0\""
"	},"
"	\"Room4.location\":{"
"		\"0\":\"_id0\""
"	},"
"	\"Optimization.location\":{"
"		\"0\":\"_id1\""
"	},"
"	\"Controller.location\":{"
"		\"0\":\"_id2\","
"		\"1\":\"Wait\","
"		\"2\":\"_id4\""
"	}"
"},\"regressors\":{"
"	\"(195)\":"
"		{\"type\":\"act->point->val\",\"representation\":\"simpletree\",\"minimize\":1,\"regressor\":"
"			{"
"			\"11\" : 1.767005809675129,"
"			\"12\" : 0.5472758743181302"
"			}"
"		},"
"	\"(180)\":"
"		{\"type\":\"act->point->val\",\"representation\":\"simpletree\",\"minimize\":1,\"regressor\":"
"			{"
"			\"11\" : 2.107357559876886,"
"			\"12\" : 0.841995896951961"
"			}"
"		},"
"	\"(225)\":"
"		{\"type\":\"act->point->val\",\"representation\":\"simpletree\",\"minimize\":1,\"regressor\":"
"			{"
"			\"11\" : 0.3895094003837855,"
"			\"12\" : 1.438914715299473"
"			}"
"		},"
"	\"(120)\":"
"		{\"type\":\"act->point->val\",\"representation\":\"simpletree\",\"minimize\":1,\"regressor\":"
"			{"
"			\"11\" : 3.562263183988571,"
"			\"12\" : 3.01531941566271"
"			}"
"		},"
"	\"(90)\":"
"		{\"type\":\"act->point->val\",\"representation\":\"simpletree\",\"minimize\":1,\"regressor\":"
"			{"
"			\"11\" : 4.595961684754502,"
"			\"12\" : 4.520572775155376"
"			}"
"		},"
"	\"(75)\":"
"		{\"type\":\"act->point->val\",\"representation\":\"simpletree\",\"minimize\":1,\"regressor\":"
"			{"
"			\"11\" : 5.200779712113728,"
"			\"12\" : 1.1120904838198"
"			}"
"		},"
"	\"(165)\":"
"		{\"type\":\"act->point->val\",\"representation\":\"simpletree\",\"minimize\":1,\"regressor\":"
"			{"
"			\"11\" : 2.422010195182813,"
"			\"12\" : 1.850311355224428"
"			}"
"		},"
"	\"(240)\":"
"		{\"type\":\"act->point->val\",\"representation\":\"simpletree\",\"minimize\":1,\"regressor\":"
"			{"
"			\"11\" : 0.1842626478060555,"
"			\"12\" : 0.5429048553483474"
"			}"
"		},"
"	\"(150)\":"
"		{\"type\":\"act->point->val\",\"representation\":\"simpletree\",\"minimize\":1,\"regressor\":"
"			{"
"			\"11\" : 1.802785410011689,"
"			\"12\" : 2.744820174794056"
"			}"
"		},"
"	\"(135)\":"
"		{\"type\":\"act->point->val\",\"representation\":\"simpletree\",\"minimize\":1,\"regressor\":"
"			{"
"			\"11\" : 4.082210512180363,"
"			\"12\" : 1.046474095029413"
"			}"
"			},"
"		\"(105)\":"
"			{\"type\":\"act->point->val\",\"representation\":\"simpletree\",\"minimize\":1,\"regressor\":"
"				{"
"				\"11\" : 3.506084076677214,"
"				\"12\" : 4.529894227895771"
"				}"
"			},"
"		\"(45)\":"
"			{\"type\":\"act->point->val\",\"representation\":\"simpletree\",\"minimize\":1,\"regressor\":"
"				{"
"				\"11\" : 2.923974662511487,"
"				\"12\" : 3.712782046133843"
"				}"
"			},"
"		\"(30)\":"
"			{\"type\":\"act->point->val\",\"representation\":\"simpletree\",\"minimize\":1,\"regressor\":"
"				{"
"				\"11\" : 1.294636921229531,"
"				\"12\" : 3.091107032150207"
"				}"
"			},"
"		\"(210)\":"
"			{\"type\":\"act->point->val\",\"representation\":\"simpletree\",\"minimize\":1,\"regressor\":"
"				{"
"				\"11\" : 1.009248500699941,"
"				\"12\" : 1.028772620412916"
"				}"
"			},"
"		\"(15)\":"
"			{\"type\":\"act->point->val\",\"representation\":\"simpletree\",\"minimize\":1,\"regressor\":"
"				{"
"				\"11\" : 2.855407148131672,"
"				\"12\" : 1.352930166302546"
"				}"
"			}"
"	}"
"}";



BOOST_AUTO_TEST_CASE(SimpleUnorderedKeyLoad)
{
    std::stringstream ss(simple_unordered_strategy);
    auto strategy = SimpleTree::parse(ss, false, false, 0);
    double disc[1] = {15.0};
    BOOST_CHECK_EQUAL(strategy.value(disc, nullptr, 11), 2.855407148131672);
    BOOST_CHECK_EQUAL(strategy.value(disc, nullptr, 12), 1.352930166302546);

    disc[0] = 210;
    BOOST_CHECK_EQUAL(strategy.value(disc, nullptr, 11), 1.009248500699941);
    BOOST_CHECK_EQUAL(strategy.value(disc, nullptr, 12), 1.028772620412916);

    disc[0] = 30;
    BOOST_CHECK_EQUAL(strategy.value(disc, nullptr, 11), 1.294636921229531);
    BOOST_CHECK_EQUAL(strategy.value(disc, nullptr, 12), 3.091107032150207);
}

BOOST_AUTO_TEST_CASE(UnorderedKeyLoad)
{
    std::stringstream ss(unordered_strategy);
    auto strategy = SimpleTree::parse(ss, false, false, 0);
    double disc[1] = {15.0};
    BOOST_CHECK_EQUAL(strategy.value(disc, nullptr, 11), 2.855407148131672);
    BOOST_CHECK_EQUAL(strategy.value(disc, nullptr, 12), 1.352930166302546);

    disc[0] = 210;
    BOOST_CHECK_EQUAL(strategy.value(disc, nullptr, 11), 1.009248500699941);
    BOOST_CHECK_EQUAL(strategy.value(disc, nullptr, 12), 1.028772620412916);

    disc[0] = 30;
    BOOST_CHECK_EQUAL(strategy.value(disc, nullptr, 11), 1.294636921229531);
    BOOST_CHECK_EQUAL(strategy.value(disc, nullptr, 12), 3.091107032150207);
}
