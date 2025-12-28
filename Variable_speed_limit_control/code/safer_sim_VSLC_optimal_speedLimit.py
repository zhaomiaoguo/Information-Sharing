# code written by Fatima Afifah
# Implementation of VSLC to improve safety and mobility


# import traci related packages
from __future__ import absolute_import
from __future__ import print_function
import os
import sys
import optparse
import random
import xml.etree.ElementTree as ET

# we need to import python modules from the $SUMO_HOME/tools directory
if 'SUMO_HOME' in os.environ:
    tools = os.path.join(os.environ['SUMO_HOME'], 'tools')
    sys.path.append(tools)
else:
    sys.exit("please declare environment variable 'SUMO_HOME'")

from sumolib import checkBinary  # noqa
from sumolib.net import readNet  # noqa
import traci  # noqa

# import learning related packages
#%matplotlib notebook

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import matplotlib.cm as cm
import pickle  # to save/load Q-Tables
import time  # using this to keep track of our saved Q-Tables.
from tqdm import notebook as tqdm
from tensorflow.keras import backend
from tensorflow.keras.models import Sequential, load_model
from tensorflow.keras.layers import Dense, Dropout, Conv2D, MaxPooling2D, Activation, Flatten
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.callbacks import TensorBoard
import tensorflow as tf
from collections import deque
from PIL import Image
#import cv2

from textwrap import wrap
import pandas as pd

# global parameters
DISCOUNT = 0.95
REPLAY_MEMORY_SIZE = 5_000  # How many last steps to keep for model training
MIN_REPLAY_MEMORY_SIZE = 500  # Minimum number of steps in a memory to start training
MINIBATCH_SIZE = 64  # How many steps (samples) to use for training
UPDATE_TARGET_EVERY = 5  # Terminal states (end of episodes)
MODEL_NAME = '2x10'
MIN_REWARD = -4750.00  # For model save
MEMORY_FRACTION = 0.20
LEARNING_RATE = 1
EFFECT_PERIODS_START = 1
EFFECT_PERIODS_END = 2

# Environment settings
EPISODES = 800

# Exploration settings
epsilon = 0.90  # not a constant, going to be decayed
EPSILON_DECAY = 0.995
MIN_EPSILON = 0.01

#  Stats settings
AGGREGATE_STATS_EVERY = 5  # episodes
SHOW_PREVIEW = False
MEMORY_CAPACITY = 64

#EP_MAX = 600
LR_A = 0.0002    # learning rate for actor
LR_C = 0.0005   # learning rate for critic
GAMMA = 0.9      # reward discount
TAU = 0.005
#TAU = 0.999# soft replacement
#MEMORY_CAPACITY = 64
MEMORY_CAPACITY = 1000000
BATCH_SIZE = 256

RENDER = False

LOAD_MODEL = None
DATA_DIR = "data_segments/"


class SUMOEnv:    
    RETURN_IMAGES = False
    ACCIDENT_START =20
    ACCIDENT_END = 80 
    ACCIDENT_EDGE_ID = "e5.66"
    INFO_EDGE_ID = "e1"
    episode_step = 0
    speed_value = 17.8816  
    speed_value1 = 2.2352
    alpha = 0
    #17.8816
    #29.06
    #24.5872
    

    max_speed = {"e0":speed_value, "e1":speed_value,"e2":speed_value,"e3":speed_value,"e4":speed_value,"e6":speed_value}
    max_speed_accident = {(ACCIDENT_EDGE_ID, 1): 0.001, (ACCIDENT_EDGE_ID, 0): 20}
    estimate_tt = {ACCIDENT_EDGE_ID: 10, "e5": 10, "e5.33": 10} # if there is no information, people will assume it takes 10 secons to travel trhough e5, e5.33, and e5.66
    
    edges = ['e0 e1 e2 e6',\
             'e0 e3 e4 e6']
    control_section = ['e1']
    state_edges = ['e0','e1','e2','e3','e4','e6']
    VSLlist = ['e1']
    EPSILON = 0.1 

     # get edge information once at the beginning of the simulation. 
    def get_edge_info():
        def get_edge_length(edge_id_list):
            edge_length = {}
            for edge_id in edge_id_list:
                edge_length[edge_id] = traci.lane.getLength(edge_id + "_0")
            return edge_length

        def get_edge_lanes(edge_id_list):
            edge_lanes = {}
            for edge_id in edge_id_list:
                edge_lanes[edge_id] = traci.edge.getLaneNumber(edge_id)
            return edge_lanes
        traci.start(["sumo", "-c", DATA_DIR+"slimer2.sumocfg"])
        # get the length of each edge, do not include the intersection connection edge
        edge_id_list = list(traci.edge.getIDList())
        edge_id_list = [x for x in edge_id_list if ":" not in x] # remove all the intersections
        edge_length = get_edge_length(edge_id_list) 
        edge_lanes = get_edge_lanes(edge_id_list)
        
        traci.close(["SUMO", "-c", DATA_DIR+"slimer2.sumocfg"])

        return edge_id_list, edge_length, edge_lanes
    
    
    #edge_id_list, edge_length, edge_lanes,control_section,state_edges,VSLlist = get_edge_info()
    edge_id_list, edge_length, edge_lanes = get_edge_info()
    
    
            
    
    
    
    
    

    OBSERVATION_SPACE_VALUES = len(state_edges) * 2  # 3: 1. ACCIDENT? 2. time step. This is finite horizon  
    ACTION_SPACE_SIZE = len(VSLlist) # share and not share information
    
    # the dict! (colors) 1,2,3 are state index
    d = {1: (255, 175, 0),
         2: (0, 255, 0),
         3: (0, 0, 255)}
    
    def set_max_speed(self, max_speed, v):
        for edge_id in self.edge_id_list: 
            traci.edge.setMaxSpeed(edgeID=edge_id, speed=max_speed[edge_id])
        number_of_edges = len(self.VSLlist)
        for j in range(number_of_edges):
            traci.edge.setMaxSpeed(self.VSLlist[j], v[j])
                
        
            

    def update_max_speed(self, edge_id, accident):
        # set the speed limit of edge e5
        traci.edge.setMaxSpeed(edgeID=edge_id, speed=self.max_speed_accident[edge_id, accident])

    

    def is_accident(self):
        if self.episode_step >= self.ACCIDENT_START and self.episode_step<self.ACCIDENT_END:
            accident = 1
        else:
            accident = 0
        
        return accident

    def get_tt_set_adapt_tt(self):
        travel_time = {}
    # current travel time for each edge
        for edge_id in self.edge_id_list:
            travel_time[edge_id] = traci.edge.getTraveltime(edge_id) # real travel time
            traci.edge.adaptTraveltime(edgeID=edge_id, time=travel_time[edge_id]) # all the edges use the real travel time except for edge e5
        
        for edge_id in self.estimate_tt.keys():
            traci.edge.adaptTraveltime(edgeID=edge_id, time=self.estimate_tt[edge_id])# people are uncertain about the travel time on e5 and assume it is 5 

        return travel_time

    def reroute_veh(self):
        veh = {}
        
        for edge_id in self.edge_id_list:
            veh[edge_id] = traci.edge.getLastStepVehicleIDs(edgeID=edge_id)
            for veh_id in veh[edge_id]:
                traci.vehicle.rerouteTraveltime(vehID=veh_id, currentTravelTimes=True)
                traci.vehicle.setRoutingMode(vehID=veh_id, routingMode= 0)
    
    
    
            
    

    # considering only speed+density as state
    def get_state(self):
        state = []
        num_vehicle = {}
        density = {}
        speed = {}
        flow = {}
        for edge_id in self.state_edges:
            num_vehicle[edge_id] = traci.edge.getLastStepVehicleNumber(edgeID=edge_id)
            density[edge_id] = num_vehicle[edge_id] / self.edge_length[edge_id]/self.edge_lanes[edge_id]
            #getLastStepVehicleIDs, if no vehicle, speed = max
            speed[edge_id] = traci.edge.getLastStepMeanSpeed(edgeID = edge_id)
            #state = state+[density[edge_id],speed[edge_id]]
            flow[edge_id] = density[edge_id]*speed[edge_id]
            state= state+[density[edge_id]]+[speed[edge_id]]
        return state



        
    
    #####################  set speed limit  #################### 
    def set_vsl(self, v):
        number_of_edges = len(self.VSLlist)
        for j in range(number_of_edges):
            traci.edge.setMaxSpeed(self.VSLlist[j], v[j][0])
            
    

    
    def calc_reward(self):    

        vidlist = traci.edge.getIDList()
        ttc = 0
        veh = {}
        for edge_id in self.edge_id_list:
            veh[edge_id] = traci.edge.getLastStepVehicleIDs(edgeID=edge_id)
            for v in veh[edge_id]:
                leader = traci.vehicle.getLeader(vehID = v,dist = 0) # just focus on its current lane
                if leader != None:
                    lead_speed = traci.vehicle.getSpeed(vehID = leader[0])
                    follow_speed = traci.vehicle.getSpeed(vehID = v)
                    speed_diff = (follow_speed - lead_speed)
                    space_gap = leader[1]+traci.vehicle.getMinGap(vehID = v) # leader has two elements: [0] vehID; [1] distance - minGap
                    if speed_diff >= self.EPSILON:
                        ttc= ttc + space_gap/speed_diff
                        
        num_vehicle = {}
        travel_time = {}
        vidlist = traci.edge.getIDList()
        total_tt = 0
        tt = {}
        for edge_id in self.edge_id_list:
            num_vehicle[edge_id] = traci.edge.getLastStepVehicleNumber(edgeID=edge_id)
            travel_time[edge_id] = traci.edge.getTraveltime(edge_id) # real travel time
            traci.edge.adaptTraveltime(edgeID=edge_id, time=travel_time[edge_id]) 
            total_tt = total_tt + travel_time[edge_id] * num_vehicle[edge_id]
             
            
        
        reward = ttc * self.alpha - (1-self.alpha) * total_tt

        return reward


    
    def reset(self, action):
        self.episode_step = 0
        traci.load(['-c', DATA_DIR+"slimer2.sumocfg"])
        self.set_max_speed(self.max_speed, action)
        
        
        traci.simulationStep() # warm up: simulate0s to set up the vehicle locations and estimate speed
        state_observe= self.get_state()
        
        return state_observe

    def step(self, action):
        self.set_max_speed(self.max_speed, action)
        self.reroute_veh()
        traci.simulationStep() # simulate one step after decision is made to get new state
        self.episode_step += 1
        new_observation= self.get_state()
        reward = self.calc_reward()
        
        done = False
        if traci.simulation.getMinExpectedNumber() == 0:
            done = True
                  
        return new_observation,reward, done, self.episode_step



    
    
env = SUMOEnv()


class TD3Agent(object):
    def __init__(self, a_dim, s_dim,):
        #self.memory = Memory(capacity=MEMORY_CAPACITY)
        self.memory = np.zeros((MEMORY_CAPACITY, s_dim * 2 + a_dim + 1), dtype=np.float32)
        self.pointer = 0
        #self.sess = tf.Session()
        self.sess = tf.compat.v1.Session()
        tf.compat.v1.disable_eager_execution()
        #tf.compat.v1.variable_scope('Actor')

        self.a_dim, self.s_dim = a_dim, s_dim
        self.S = tf.compat.v1.placeholder(tf.float32, [None, s_dim], 's')
        self.S_ = tf.compat.v1.placeholder(tf.float32, [None, s_dim], 's_')
        self.R = tf.compat.v1.placeholder(tf.float32, [None, 1], 'r')

        self.a = self._build_a(self.S,)  
        q = self._build_c(self.S, self.a, )
        a_params = tf.compat.v1.get_collection(tf.compat.v1.GraphKeys.TRAINABLE_VARIABLES, scope='Actor')
        c_params = tf.compat.v1.get_collection(tf.compat.v1.GraphKeys.TRAINABLE_VARIABLES, scope='Critic')
        ema = tf.train.ExponentialMovingAverage(decay=1 - TAU)          # soft replacement

        def ema_getter(getter, name, *args, **kwargs):
            return ema.average(getter(name, *args, **kwargs))

        self.update_targetActor = ema.apply(a_params)
        self.update_targetCritic = ema.apply(c_params)
        target_update = [self.update_targetActor, self.update_targetCritic]
        #target_update = [ema.apply(a_params), ema.apply(c_params)]      # soft update operation
        a_ = self._build_a(self.S_, reuse=True, custom_getter=ema_getter)   # target network for actor
        q_ = self._build_c(self.S_, a_, reuse=True, custom_getter=ema_getter)  # target network for critic

        a_loss = - tf.reduce_mean(q)  # maximize the q
        self.atrain = tf.compat.v1.train.AdamOptimizer(LR_A).minimize(a_loss, var_list=a_params)
        self.td = self.R + GAMMA * q_ - q

        with tf.control_dependencies(target_update):    # soft replacement happened at here
            q_target = self.R + GAMMA * q_ 
            td_error = tf.compat.v1.losses.mean_squared_error(labels=q_target, predictions=q)
            self.ctrain = tf.compat.v1.train.AdamOptimizer(LR_C).minimize(td_error, var_list=c_params)

        self.sess.run(tf.compat.v1.global_variables_initializer())
        self.saver = tf.compat.v1.train.Saver(max_to_keep = 1)
        
    
    def choose_action(self, s):
        return self.sess.run(self.a, {self.S: s[np.newaxis, :]})[0]
    
    
    def learn(self):
        record_range = min(self.pointer, MEMORY_CAPACITY)
        indices = np.random.choice(record_range, size=BATCH_SIZE)
        bt = self.memory[indices, :]
        bs = bt[:, :self.s_dim]
        ba = bt[:, self.s_dim: self.s_dim + self.a_dim]
        br = bt[:, -self.s_dim - 1: -self.s_dim]
        bs_ = bt[:, -self.s_dim:]

        self.sess.run([self.atrain, self.update_targetActor], {self.S: bs})
        self.sess.run([self.ctrain, self.update_targetCritic], {self.S: bs, self.a: ba, self.R: br, self.S_: bs_})


#    def store_transition(self, s, a, r, s_):
#        transition = np.hstack((s, a, r, s_))
#        self.memory.store(transition)
    # replay buffer
    def store_transition(self, s, a, r, s_):
        transition = np.hstack((s, a, [r], s_))
        index = self.pointer % MEMORY_CAPACITY  # replace the old memory with new memory
        self.memory[index, :] = transition
        self.pointer += 1


    def _build_a(self, s, reuse=None, custom_getter=None):
        trainable = True if reuse is None else False
        with tf.compat.v1.variable_scope('Actor', reuse=reuse, custom_getter=custom_getter):
            neta = tf.compat.v1.layers.dense(s, 60, activation=tf.nn.relu, name='l1', trainable=trainable)
            a = tf.compat.v1.layers.dense(neta, self.a_dim, activation=tf.nn.sigmoid, name='l2', trainable=trainable,  use_bias=False)
            return tf.multiply(a, 8, name='scaled_a')

    def _build_c(self, s, a, reuse=None, custom_getter=None):
        trainable = True if reuse is None else False
        with tf.compat.v1.variable_scope('Critic', reuse=reuse, custom_getter=custom_getter):
            n_l1 = 50
            w1_s = tf.compat.v1.get_variable('w1_s', [self.s_dim, n_l1], trainable=trainable)
            w1_a = tf.compat.v1.get_variable('w1_a', [self.a_dim, n_l1], trainable=trainable)
            b1 = tf.compat.v1.get_variable('b1', [1, n_l1], trainable=trainable)
            netc = tf.nn.relu(tf.matmul(s, w1_s) + tf.matmul(a, w1_a) + b1)
            return tf.compat.v1.layers.dense(netc, 1, trainable=trainable)
    
    def savemodel(self,):
        self.saver.save(self.sess,'ddpg_networkss_withoutexplore/' + 'ddpg.ckpt')
        
    def loadmodel(self,):
        loader = tf.train.import_meta_graph('ddpg_networkss_withoutexplore/ddpg.ckpt.meta')
        loader.restore(self.sess, tf.train.latest_checkpoint("ddpg_networkss_withoutexplore/"))
env = SUMOEnv()
agent = TD3Agent(s_dim = env.OBSERVATION_SPACE_VALUES, a_dim = env.ACTION_SPACE_SIZE)


def from_a_to_mlv(a):
    return 13.4112 + 2.2352*np.floor(a)

# To store reward history of each episode
ep_reward_list = []
# To store average reward history of last few episodes
avg_reward_list = []
avg_reward_list20 = []


all_ep_r=[]
all_ttc =[]
all_total_tt =[]
aggr_ep_rewards = {'ep': [], 'avg': [], 'max': [], 'min': []}
aggr_ep_ttc = {'ep': [], 'avg': [], 'max': [], 'min': []}
aggr_ep_total_ttc = {'ep': [], 'avg': [], 'max': [], 'min': []}
VSL_actn = {}
vsl_r={}

#VSL_actn = {'ep': [], 'step': [], 'edge': [], 'speed_limit': []}
# For more repetitive results
random.seed(3)
np.random.seed(4)
tf.random.set_seed(5)



# Create models folder
if not os.path.isdir('models'):
    os.makedirs('models')

# start sumo first 
traci.start(["sumo", "-c", DATA_DIR+"slimer2.sumocfg"])
ttc = 0
total_tt = 0
total_step = 0
stime = np.zeros(env.OBSERVATION_SPACE_VALUES,)




# Iterate over episodes
for episode in tqdm.tqdm(range(1, EPISODES + 1), ascii=True, unit='episodes'):
    
    state_history = {}
    accident_history = {}
    action_history = {}
    reward_history = {}
    
    

    
    ep_r = 0
    #episode_tt = 0
    step = 0
    #VSL_actn[episode, step] = []
    action = env.speed_value*np.ones(env.ACTION_SPACE_SIZE,)
    # Reset environment and get initial state
    current_state = env.reset(action) # current state are all the density and speed of each edge (including intersection?)
    state_history[step] = current_state
    stime[0:env.OBSERVATION_SPACE_VALUES] = current_state

    var = 1.5
    
    
    
    done = False
    while not done:
        action = agent.choose_action(stime)
        action= np.clip(np.random.laplace(action, var), 0, 7.99) # 40 to 75mph
        
        
        
        v = from_a_to_mlv(action)
        stime_ = np.zeros(env.OBSERVATION_SPACE_VALUES,)
        print(v)
        
        number_of_edges = len(env.VSLlist)
        
        
        for j in range(number_of_edges):
            VSL_actn[episode, step, env.VSLlist[j]] = v[j]
            #print(v[j])
            
        s_,r, done ,simulationSteps= env.step(v)
        
        

        
        stime_[0:env.OBSERVATION_SPACE_VALUES] = s_
        #stime_[env.OBSERVATION_SPACE_VALUES-1] = simulationSteps/18000
        agent.store_transition(stime, action, r, stime_)
        total_step = total_step + 1

        if total_step > MEMORY_CAPACITY:
            agent.learn()
        stime = stime_
        ep_r += r
        vsl_r[episode, step]=r
        step += 1
        #print("ep_r", ep_r)
    all_ep_r.append(ep_r)
    #print("all_ep_r",all_ep_r)
    if not episode % AGGREGATE_STATS_EVERY or episode == 1:
        average_reward = sum(all_ep_r[-AGGREGATE_STATS_EVERY:])/len(all_ep_r[-AGGREGATE_STATS_EVERY:])
        min_reward = min(all_ep_r[-AGGREGATE_STATS_EVERY:])
        max_reward = max(all_ep_r[-AGGREGATE_STATS_EVERY:])
        aggr_ep_rewards['ep'].append(episode)
        aggr_ep_rewards['avg'].append(average_reward)
        aggr_ep_rewards['max'].append(max_reward)
        aggr_ep_rewards['min'].append(min_reward)
        
    ep_reward_list.append(ep_r)

    # Mean of last 40 episodes
    avg_reward = np.mean(ep_reward_list[-40:])
    avg_reward2 = np.mean(ep_reward_list[-20:])
    print("Episode * {} * Avg Reward is ==> {}".format(episode, avg_reward))
    avg_reward_list.append(avg_reward)
    avg_reward_list20.append(avg_reward2)
    
    
        
traci.close(["SUMO", "-c", DATA_DIR+"slimer2.sumocfg"])    
agent.savemodel()


rewards = []
for e, s, ed, sp in zip(aggr_ep_rewards['ep'],aggr_ep_rewards['avg'],aggr_ep_rewards['max'],aggr_ep_rewards['min']):
    data2 = { 'episode': e, 'avg': s, 'max': ed, 'min': sp }
    rewards.append(data2)

reward = pd.DataFrame(rewards)

reward.to_csv('reward_exp_lr_epgr_up.csv')


ep =[]
step =[]
edge = []
speed_limit=[]
for k, v in VSL_actn.items():
    ep.append(k[0])
    step.append(k[1])
    edge.append(k[2])
    speed_limit.append(v)

data = []
for e, s, ed, sp in zip(ep,step,edge,speed_limit):
    data1 = { 'episode': e, 'step': s, 'edge': ed, 'speed_limit': sp }
    data.append(data1)

df = pd.DataFrame(data)
df['speed_limit_mph'] = df['speed_limit'] * 2.23694
df.to_csv('action_exp_lr_epgr_up.csv')

ep = list(range(1,EPISODES+1))
dict = {"Episode": ep, "reward": avg_reward_list}
df2 = pd.DataFrame(dict)
df2.to_csv('avg_reward_exp_lr_epgr_up.csv')
dict2 = {"Episode": ep, "reward": avg_reward_list20}
df3 = pd.DataFrame(dict2)
df3.to_csv('avg_reward20_exp_lr_epgr_up.csv')
                                           
