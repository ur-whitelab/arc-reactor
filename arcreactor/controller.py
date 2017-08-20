import zmq
import zmq.asyncio
import time
import argparse
import asyncio
from .server import start_server


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
        print('listening for topic: {}'.format(sub_topic))
        self.vision_sock.subscribe(sub_topic.encode())

        #register publishing socket
        zmq_uri = 'tcp://{}:{}'.format(cc_hostname, zmq_pub_port)
        print('Connecting PUB Socket to {}'.format(zmq_uri))
        self.pub_sock = self.ctx.socket(zmq.PUB)
        self.pub_sock.connect(zmq_uri)

        self.frequency = 1
        self.simulation_state = 'test'

    async def handle_start(self,server_port):
        '''Begin processing reactor simulation'''

        start_server(self,server_port)
        print('Started arcreactor server')
        import sys
        sys.stdout.flush()

        while True:
            await self.update_loop()

    async def update_simulation(self, vision_state):
        asyncio.sleep(0.25)
        pass

    async def update_loop(self):
        msg_parts = await self.vision_sock.recv_multipart()
        vision_update = msg_parts[1]
        start_time = time.time()
        await self.update_simulation(vision_update)
        #exponential moving average of update frequency
        self.frequency = self.frequency * 0.8 +  0.2 / max(0.0001, time.time() - start_time)
        await self.pub_sock.send_multipart(['simulation-update'.encode(), self.simulation_state.encode()])


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