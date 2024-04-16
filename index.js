const express = require('express')
const app = express()

app.get('/api/hello/:name', (req, res)=>{
    var name =req.params.name
    res.send(
    `<!DOCTYPE html>
    <html>
    <head>
      <title>B1TMET Message</title>
      <link rel="icon" type="image/x-icon" href="/images/favicon.png">
    </head>
    <body>
    
    <!-- <h1>This is a Heading</h1>
    <p>This is a paragraph.</p> -->
    
    <h1 style="text-align: center;"><strong>B1TMET HTML Test Message</strong></h1>
    <h2 style="text-align: center;"><strong>&nbsp;</strong></h2>
    <h2 style="text-align: center;"><strong>Hi&nbsp;<u> ${name}</u>, this is a &nbsp;test message. &nbsp;</strong></h2>
    <h2 style="text-align: center;"><strong>&nbsp;</strong></h2>
    <h2 style="text-align: center;"><strong>Best regards, from B1TMET team ;)</strong></h2>
    
    </body>
    </html>`
    )
    console.log(`Hello ${name}!`)
})





var puerto = process.env.PORT ?? 4000
app.listen(puerto, ()=>{
    require('log-timestamp')
    console.log('#########################################################')
    console.log(`Server Client Project on port ${puerto}: Connected`);
    console.log('#########################################################\n')

});

