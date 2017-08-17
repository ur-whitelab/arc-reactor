import zmq
import zmq.asyncio
import time
import argparse
import asyncio
from .server import start_server


zmq.asyncio.install()

class Controller:
    '''Controls flow of reactor program'''
    def __init__(self, zmq_uri):
        self.ctx = zmq.asyncio.Context()
        print('Connecting SUB Socket to {}'.format(zmq_uri))
        self.vision_sock = self.ctx.socket(zmq.SUB)
        self.vision_sock.connect(zmq_uri)
        sub_topic = 'vision-update'
        print('listening for topic: {}'.format(sub_topic))
        self.vision_sock.subscribe(sub_topic)
        self.frequency = 1

    async def handle_start(self,server_port):
        '''Begin processing reactor simulation'''

        start_server(self,server_port)
        print('Started arcreactor server')
        import sys
        sys.stdout.flush()

        #need to discover vision publisher information        

        while True:
            await self.update_loop()

    def update_simulation(self):
        pass

    async def update_loop(self):        
        vision_state = await self.vision_sock.recv_multipart()
        start_time = time.time(vision_state)
        self.update_simulation()
        #exponential moving average of update frequency
        self.frequency = self.frequency * 0.8 +  0.2 / (time.time() - start_time)


def init(video_filename, server_port, zmq_port, hostname):
    c = Controller('tcp://{}:{}'.format(hostname, zmq_port))
    asyncio.ensure_future(c.handle_start(server_port))
    loop = asyncio.get_event_loop()
    loop.run_forever()


def main():
    parser = argparse.ArgumentParser(description='Project ARC ')
    parser.add_argument('--server-port', help='port to run server', default='8888', dest='server_port')
    parser.add_argument('--zmq-port', help='port for publishing zmq update', default=5000, dest='zmq_port')
    parser.add_argument('--hostname', help='hostname for process', default='*')
    args = parser.parse_args()
    init(args.server_port, args.zmq_port, args.hostname)