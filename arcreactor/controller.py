import zmq
import zmq.asyncio
import sys
import time
import argparse
import asyncio
from server import start_server
from protobufs.graph_pb2 import Graph
from protobufs.kinetics_pb2 import SystemKinetics
from analysis import Analyzer
from simulation import Simulation

import copy
zmq.asyncio.install()

class Controller:
    '''Controls flow of reactor program'''
    def __init__(self, zmq_sub_port, zmq_pub_port, cc_hostname):
        self.ctx = zmq.asyncio.Context()

        #subscribe to publishing socket
        zmq_uri = 'tcp://{}:{}'.format(cc_hostname, zmq_sub_port)
        print('Connecting SUB Socket to {}'.format(zmq_uri))
        self.vision_sock = self.ctx.socket(zmq.SUB)
        self.vision_sock.connect(zmq_uri)
        #we only want vision updates
        sub_topic = 'vision-update'
        print('listening for topic: {} on SUB'.format(sub_topic))
        self.vision_sock.subscribe(sub_topic)

        #register publishing socket
        zmq_uri = 'tcp://{}:{}'.format(cc_hostname, zmq_pub_port)
        print('Connecting PUB Socket to {}. Will publish simulation-update'.format(zmq_uri))
        self.pub_sock = self.ctx.socket(zmq.PUB)
        self.pub_sock.connect(zmq_uri)

        self.analyzer = Analyzer()
        self.frequency = 1
        self.stream_names = self.analyzer.stream_names

        self.graph = Graph()
        self.simulator = Simulation(0)
        self.simulation_state = SystemKinetics()

    async def handle_start(self,server_port):
        '''Begin processing reactor simulation'''

        start_server(self,server_port)
        print('Started arcreactor server')
        sys.stdout.flush()

        while True:
            await self.update_loop()
            sys.stdout.flush()

    async def update_simulation(self):
        '''Update simulation every time a message is received from vision'''

        if self.simulator.start_time == 0:
            self.simulator.start_time = self.graph.time
        self.simulation_state.time = self.graph.time
        new_graph = copy.copy(self.graph)
        new_sim_state = copy.copy(self.simulation_state)
        self.simulation_state = await self.simulator.calculate(new_sim_state, new_graph)
        #print('Called calculate() in update_simulation(). Now self.simulation_state is {}'.format(self.simulation_state))
        #print('and self.graph was {}'.format(self.graph))
        await asyncio.sleep(0)
        return self.simulation_state

    async def update_loop(self):
        '''Send messages to visualization'''
        msg_parts = await self.vision_sock.recv_multipart()
        self.graph.ParseFromString(msg_parts[1])
        start_time = time.time()
        await self.update_simulation()
        #exponential moving average of update frequency
        self.frequency = self.frequency * 0.8 +  0.2 / max(0.0001, time.time() - start_time)
        await self.pub_sock.send_multipart(['simulation-update'.encode(), self.simulation_state.SerializeToString()])


def init(server_port, zmq_sub_port, zmq_pub_port, cc_hostname):
    c = Controller(zmq_sub_port, zmq_pub_port, cc_hostname)
    asyncio.ensure_future(c.handle_start(server_port))
    loop = asyncio.get_event_loop()
    loop.run_forever()


def main():
    parser = argparse.ArgumentParser(description='Project ARC ')
    parser.add_argument('--server-port', help='port to run server', default='8888', dest='server_port')
    parser.add_argument('--zmq-sub-port', help='port for receiving zmq sub update', default=5000, dest='zmq_sub_port')
    parser.add_argument('--cc-hostname', help='hostname for cc to receive zmq pub updates', default='localhost', dest='cc_hostname')
    parser.add_argument('--zmq-pub-port', help='port for publishing my zmq updates', default=2400, dest='zmq_pub_port')
    args = parser.parse_args()
    init(args.server_port, args.zmq_sub_port, args.zmq_pub_port, args.cc_hostname)