
#!/usr/bin/env python3
'''Minimal sevrer to send/receive updates via HTTP instead of ZMQ'''

import tornado.web
from tornado.platform.asyncio import AsyncIOMainLoop
import asyncio
import os
import json

AsyncIOMainLoop().install()

RESOURCES = os.path.join(os.path.dirname(__file__), os.pardir, 'resources')
WEB_STRIDE = 1

class HtmlPageHandler(tornado.web.RequestHandler):
    async def get(self, file_name='index.html'):
        # Check if page exists
        www = os.path.join(RESOURCES, file_name)
        if os.path.exists(www):
            # Render it
            self.render(www)
        else:
            # Page not found, generate template
            err_tmpl = tornado.template.Template("<html> Err 404, Page {{ name }} not found</html>")
            err_html = err_tmpl.generate(name=file_name)
            # Send response
            self.finish(err_html)


class StatsHandler(tornado.web.RequestHandler):
    '''Provides stats on speed of processing'''
    def initialize(self, controller):
        self.controller = controller

    async def get(self):
        self.set_header("Access-Control-Allow-Origin", "*")
        self.set_header("Access-Control-Allow-Headers", "x-requested-with")
        self.set_header('Access-Control-Allow-Methods', 'GET, OPTIONS')
        self.write(json.dumps(self.controller.__dict__, default=lambda x: ''))


def start_server(controller, port=8888):
    app = tornado.web.Application([
        (r"/",HtmlPageHandler),
        (r"/stats", StatsHandler, {'controller': controller})
    ])
    print('Starting server on port {}'.format(port))
    app.listen(port)

async def get_process_info(process, url):
    http_client = AsyncHTTPClient()
    response = await http_client.fetch(hostname + ':' + port + '/process/' + process)
    return json.loads(response.body)
