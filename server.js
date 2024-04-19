const express = require('express')
const mod_fs = require('fs')
const app = express()
const puerto = process.env.PORT ?? 4000
const path = require('path')
const { exec } = require('child_process');

const bodyParser = require('body-parser');
const bodyParserErrorHandler = require('express-body-parser-error-handler')


// const fecha = new Date();



// //// Log a Archivo 1

// var fs = require('fs'); var util = require('util');
// var log_file = fs.createWriteStream(__dirname + '/server_CLientProject.log', {flags : 'a'});
//                                                                                 // Or 'w' to truncate the file every time the process starts.

// var log_stdout = process.stdout;

// console.log = function(d) { //
//     log_file.write(util.format(d) + '\n');
//     log_stdout.write(util.format(d) + '\n');
// };

// ////




// app.use(express.static(__dirname+'/assets/'))
app.use(express.static(__dirname+'/assets'))

app.use(express.json({limit: "50mb"}));
app.use(express.urlencoded({ limit: "50mb", extended: true, parameterLimit: 50000 }));


app.use(bodyParser.json({limit: '50mb'}))
app.use(bodyParser.urlencoded({extended: true, limit: '50mb'}))
app.use(bodyParserErrorHandler())



app.get("/", (req, res)=>{
    // console.log(`remoteAddress: ${req.socket?.remoteAddress}`)
    const parseIp = (req) =>
    req.headers['x-forwarded-for']?.split(',').shift()
    || req.socket?.remoteAddress

    console.log(`Access detected from ${parseIp(req)}`)
    res.redirect("/api/ClientProject/")
    
})


// // Redirect function put into html code
// <script>
// var timeout = 5000;

// setTimeout(function () {
//     window.location = "/api/ClientProject/home";
//     }, timeout);
// </script>


app.get("/api/ClientProject", (req, res)=>{
    res.send(`

        <!DOCTYPE html>
        <html lang="en">

        <head>
            <meta charset="UTF-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <meta http-equiv = "refresh" content = "6; url = /api/ClientProject/home" />

            <title>Hey there!!! </title>
        </head>

        <body>
            <center>
                <form style="border:orange; border-width:5px; border-style:solid; width: 90%;">

                <h1> Hey people ;)</h1>
                <h2> In a few seconds you will be redirected to Home Page.</h2>
                <!-- <h2> Server is done coming soon</h2> -->
                <h2> Server will be ready soon </h2>

                <br>
                <h3>Best regards, <p style="color: orange;"> CFA team</p></h3>
                <br>
                <h5>See you ;)</h5>
                </form>
                    
            </center>
        </body>

        </html>

    `)

 
    // res.redirect('/api/ClientProject/home');

})



app.get("/api/ClientProject/home", (req, res)=>{
    res.sendFile("homePage.html", { root: __dirname })
    
})

app.get("/api/ClientProject/index_query_false", (req, res)=>{

    res.sendFile("indexPage.html", { root: __dirname })

})


// app.get("/api/ClientProject/index_query_true", (req, res)=>{
//     var cad = req.query.cad.replace(/['"]+/g, '')

//     let list_bool = []
//     // console.log(cad.includes('def'))
//     const cmdline = 'latiS=5.0&longW=-94.6&latiN=29.0&longE=-50.0&modelo=3&indice=1&periodoS=1&annoIRef=1961&annoIEva=1961&annoFRef=1980&annoFEva=1980&moddro=-1.0&sevdro=-1.5&extdrou=-2&umbral=-0.5&def=2&comportamiento=0&lsm=0'
//     lista_cmdline = cmdline.split('&')
//     for (let index = 0; index < lista_cmdline.length; index++) {
//         const arg = lista_cmdline[index].split('=')[0];
//         // console.log(`${arg}, ${cad.includes(arg)}`)
//         list_bool.push(cad.includes(arg+'='))
//     }
//     // console.log(list_bool)
//     // console.log('***********')

//     if(list_bool.includes(false)){
//         console.log(`Error de linea => Formato incorrecto`)

//         res.send(`
//                     <!DOCTYPE html>
//                     <html lang="en">

//                     <head>
//                         <meta charset="UTF-8">
//                         <meta name="viewport" content="width=device-width, initial-scale=1.0">
//                         <title>Ya lo tengo </title>
//                     </head>

//                     <body>
//                         <center>
//                             <h1>Error en la línea de comando!</h1>
//                             <br>
//                             <h2>Formato incorrecto</h2>
//                         </center>
//                     </body>

//                     </html>

//                     `)
//     }else{


//         const modeloi = cad.indexOf('modelo'); 
//         const modelof = modeloi+'modelo'.length
//         const modelo = cad.slice(modelof+1, modelof+2)

//         const comportamientoi = cad.indexOf('comportamiento');
//         const comportamientof = comportamientoi+'comportamiento'.length
//         const comportamiento = cad.slice(comportamientof+1, comportamientof+2)


//         if (comportamiento == '0'){
//             console.log('/api/ClientProject/query?'+cad)
//             res.redirect('/api/ClientProject/query?'+cad)

//         }else{
//             console.log(`Error de argumento => comportamiento=${comportamiento}`)
//             res.send(`
            
//             <!DOCTYPE html>
//             <html lang="en">

//             <head>
//                 <meta charset="UTF-8">
//                 <meta name="viewport" content="width=device-width, initial-scale=1.0">
//                 <title>Ya lo tengo </title>
//             </head>

//             <body>
//                 <center>
    
//                     <h1>El valor de Comportamiento debe ser cero (0) por el momento.</h1> 
                    
//                     <h2>Estoy recibiendo comportamiento=${comportamiento}</h2>
//                 </center>
//             </body>

//             </html>
//                     `)
//         }   
//     }

// })



app.get("/api/ClientProject/index_query_true", (req, res)=>{
    var cad = req.query.cad.replace(/['"]+/g, '')

    let list_bool = []
    // console.log(cad.includes('def'))
    const cmdline = 'latiS=5.0&longW=-94.6&latiN=29.0&longE=-50.0&modelo=3&indice=1&periodoS=1&annoIRef=1961&annoIEva=1961&annoFRef=1980&annoFEva=1980&moddro=-1.0&sevdro=-1.5&extdrou=-2&umbral=-0.5&def=2&comportamiento=0&lsm=0'
    lista_cmdline = cmdline.split('&')
    for (let index = 0; index < lista_cmdline.length; index++) {
        const arg = lista_cmdline[index].split('=')[0];
        // console.log(`${arg}, ${cad.includes(arg)}`)
        list_bool.push(cad.includes(arg+'='))
    }
    // console.log(list_bool)
    // console.log('***********')

    if(list_bool.includes(false)){
        console.log(`Error de linea => Formato incorrecto`)

        res.send(`
                    <!DOCTYPE html>
                    <html lang="en">

                    <head>
                        <meta charset="UTF-8">
                        <meta name="viewport" content="width=device-width, initial-scale=1.0">
                        <title>Error </title>
                    </head>

                    <body>
                        <center>
                            <h1>Error en la línea de comando!</h1>
                            <br>
                            <h2>Formato incorrecto</h2>
                        </center>
                    </body>

                    </html>

                    `)
    }else{


        if ( (cad.includes('comportamiento=0') || cad.includes('comportamiento=1'))  && cad.includes('modelo=3') && cad.includes('lsm=0')){
            console.log('/api/ClientProject/query?'+cad)
            res.redirect('/api/ClientProject/query?'+cad)

        }else{
            console.log(`Error de argumento`)
            res.send(`
            
            <!DOCTYPE html>
            <html lang="en">

            <head>
                <meta charset="UTF-8">
                <meta name="viewport" content="width=device-width, initial-scale=1.0">
                <title>Error </title>
            </head>

            <body>
                <center>
                    <h1>Error de argumentos</h1>
                    <p style="font-size: large;">El valor de <b>comportamiento</b> debe ser <b>cero (0)</b> o <b>uno (1)</b> por el momento.</p>
                    <br> 
                    <p style="font-size: large;">El valor de <b>modelo</b> debe ser <b>tres (3)</b> por el momento.</p> 
                    <br> 
                    <p style="font-size: large;">El valor de <b>lsm</b> debe ser <b>cero (0)</b> por el momento.</p> 
                    <br>
                    <br>
                    <h2> Revise estos parámetros</h2>
                    
                </center>
            </body>

            </html>
                    `)
        }   
    }

})







// app.get("/api/home", (req, res)=>{
//     res.sendFile("homePage.html")
    
// })


// app.get("/", (req, res)=>{
//     res.sendFile("./index.html")
    
// })

app.get("/api/ClientProject/query", (req, res)=>{

    const parseIp = (req) =>
    req.headers['x-forwarded-for']?.split(',').shift()
    || req.socket?.remoteAddress


    console.log(' ')
    console.log('------------------------------------------------------------')
    console.log(`IP: ${parseIp(req)}`)

    // AQUI VA EL SCRIPT DE PYTHON PARA BORRAR LOS EXISTENTES SEGUN IP
    var original_query = req.originalUrl
    var query_index = original_query.indexOf('?');
    var query_string = (query_index>=0)?req.originalUrl.slice(query_index+1):''
    query_string = JSON.stringify(query_string)
    console.log(`stdin: ${query_string}\n`)
    // res.redirect('/python_script/'+query_string);
    res.redirect('/api/ClientProject/python_script/'+query_string+'/'+parseIp(req));




    // let raw_name = req.query.name
    // let raw_fun = req.query.function

    // const name = raw_name.split(' ').join('_')
    // const fun = raw_fun.split(' ').join('_')

    // console.log(`name: ${name} fun: ${fun}`)

    // let file_down = false 
    // let data_not_found = false 

    // try {
    //     data_not_found = mod_fs.readFileSync("./notfound.html", {encoding: 'utf-8', flag: 'r'})
    // } catch (error) {
    //     res.status(200).send("Archivo [notfound.html] no encontrado. " + error)
        
    // }


    // if (name.length <= 0 || fun.length <= 0) {
    //     return res.status(200).send(data_not_found)   
    // }else{

    //     exec(`python child-processes/old_script.py ${fun} ${name}`, (error, stdout, stderr) => {

    //         if (error) {
    //         // console.log(`error: ${error.message}`);
    //         return res.status(200).send('fallos en los scripts')

    //         }
        
    //         if (stderr) {
    //         // console.log(`stderr: ${stderr}`);
    //         return res.status(200).send('fallos en los scripts')

    
    //         }
        
    //         // console.log(`stdout: ${stdout.toString()}`);

    //         try {
    //             file_down = mod_fs.readFileSync("./filedown.html", {encoding: 'utf-8', flag: 'r'})
    //             file_down = file_down.replaceAll("%namefig%", stdout)
    
    //         } catch (error) {
    //             res.status(200).send("Archivo [filedown.html] no encontrado. " + error)
                
    //         }

    //         res.status(200).send(file_down)
            

    //         });
    // }
});

app.get('/api/ClientProject/python_script/:query/:ip',(req, res)=>{
    // console.log(`aqui empieza ${req.params.query.toString()}`)
    let file_down = false 

    try{
        // const query = req.params.query;
        // const namefig = req.params.nameFig;
        
        var query = req.params.query.toString().substring(1,req.params.query.toString().length-1)
        var arg = `${req.params.query}`
        var ip = `${req.params.ip}`

        // console.log('^^^^^^^^^^^^^^^^^^^^^^^^^^^^')
        // console.log(`req.params.query: ${req.params.query}`)
        // console.log(`query: ${query}`)
        // console.log(`arg: ${arg}`)
        // console.log('^^^^^^^^^^^^^^^^^^^^^^^^^^^^')

        // console.log(`este es el argumento: ${arg}`)
        // exec(`python child-processes/script.py ${arg}`, (error, stdout, stderr) => {
        exec(`python child-processes/server_script.py ${arg} ${ip}`,{maxBuffer: 50*1024*1024}, (error, stdout, stderr) => {
            if (error) {
            console.error(`error: ${error.message}`);
            res.send(error.message)
     
            return;
            }
        
            if (stderr) {
            console.error(`stderr: ${stderr}`);
            res.send(stderr)

            return;
            }

            console.log(`stdout: ${stdout}`);

            try {
                file_down = mod_fs.readFileSync("./filedown.html", {encoding: 'utf-8', flag: 'r'})
                file_down = file_down.replaceAll("%server_name%", `${stdout}`)
                file_down = file_down.replaceAll("%client_name%", `${stdout.split('-')[1]}`)
                file_down = file_down.replaceAll("%client_ip%", `${ip}`)
                file_down = file_down.replaceAll("%client_query%", `${query}`)


            } catch (error) {
                res.status(200).send("Archivo no encontrado. " + error)
                
            }
        

            // var desired_code = '<div>Desired Code</div>'

            // mod_fs.readFileSync('./filedown1.html', 'utf8', function (err,data) {
            //   if (err) {
            //     return console.log(err);
            //   }
            
            //   //Here I need to load all the body tag content, but how?
            //   var result = data.replace('<body>html tags to be replaced</body>', desired_code);
            
            //   mod_fs.writeFileSync('./filedown1.html', result, 'utf8', function (err) {
            //      if (err) return console.log(err);
            //   });
            // });

            res.status(200).send(file_down)
            // res.status(200).send(`${output_script}`)

        })
    } catch (error) {
        res.status(500).json(error);
    }
          
} );


// app.get("/descargar/:server_name/:client_query/:client_ip/", (req, res)=>{
    app.get("/api/ClientProject/descargar/:server_name/:client_query/:client_ip", (req, res)=>{
    let file_notfound = false 

    var server_name = req.params.server_name
    var client_ip = req.params.client_ip
    var client_query = req.params.client_query
    console.log(`estoy en descargar/ y esta es la query ${client_query}`)
    var client_name = server_name.split('-')[1]
    console.log(`DOWNLOADING ${server_name} ... `)
    console.log(`client_name: ${client_name}`)

    
    res.download(__dirname+'/uploads/targz_files/'+server_name, client_name,(err)=>{
       if(err) {
           console.log(err)
       }

       exec(`rm ${__dirname}/uploads/targz_files/${server_name}`, {maxBuffer: 5 * 1024 * 1024}, async(error, stderr) => {

        if (error) {
        console.log(`rm error: ${error.message}`);


        file_notfound = mod_fs.readFileSync("./notfound.html", {encoding: 'utf-8', flag: 'r'})
        file_notfound = file_notfound.replaceAll("%client_query%", `${client_query}`)

        res.status(200).send(file_notfound)

        // res.status(200).sendFile(`${__dirname}/notfound.html`)

        // return res.status(200).send('No hay ficheros disponibles')
        // return error.toString()
        }
    
        if (stderr) {
        console.log(`rm stderr: ${stderr}`);
        // return res.status(200).send('fallos en los scripts')
        // return stderr.toString()
        }

        console.log(`${server_name} was deleted`)
        // res.redirect('/')

        // console.log('ya se borraron')
        // return stdout.toString()
        });




    });

});




app.listen(puerto, ()=>{
    require('log-timestamp')
    console.log('#########################################################')
    console.log(`Server Client Project on port ${puerto}: Connected`);
    console.log('#########################################################\n')



});


