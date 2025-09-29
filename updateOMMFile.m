function updateOMMFile(inputFile, outputFile, elements)
% updateOMMFile  Update Keplerian orbital elements in an OMM XML file
%
%   updateOMMFile(inputFile, outputFile, elements)
%
%   INPUTS:
%     inputFile   - path to original OMM XML file
%     outputFile  - path to save updated file
%     elements    - struct with fields:
%                      meanMotion [rev/day] (required)
%                      eccentricity (required)
%                      inclination [deg] (required)
%                      raan [deg] (required)
%                      argPeriapsis [deg] (required)
%                      meanAnomaly [deg] (required)
%                      epoch (optional) - datetime object or string

    if ~exist(inputFile, 'file')
        error('Input file not found: %s', inputFile);
    end

    % Validate required fields
    requiredFields = {'meanMotion', 'eccentricity', 'inclination', 'raan', 'argPeriapsis', 'meanAnomaly'};
    for i = 1:length(requiredFields)
        if ~isfield(elements, requiredFields{i}) || isempty(elements.(requiredFields{i}))
            error('%s must be provided in elements struct.', requiredFields{i});
        end
    end

    mu_earth = 3.986004418e14; % m^3/s^2
    earth_radius = 6378.137; % km
    
    % Calculate semi-major axis and derived parameters from mean motion
    n_rad_s = elements.meanMotion * 2*pi / 86400; % rev/day -> rad/s
    a_m = (mu_earth / n_rad_s^2)^(1/3);          % semi-major axis in meters
    a_km = a_m / 1000;                           % semi-major axis in km
    
    % Calculate apoapsis and periapsis altitudes
    apoapsis_alt = a_km * (1 + elements.eccentricity) - earth_radius;  % km altitude
    periapsis_alt = a_km * (1 - elements.eccentricity) - earth_radius; % km altitude
    
    % Calculate period
    period_min = 2*pi / n_rad_s / 60;            % minutes
    
    try
        xmlDoc = xmlread(inputFile);

        % Update EPOCH to match simulation time
        if isfield(elements, 'epoch') && ~isempty(elements.epoch)
            if isa(elements.epoch, 'datetime')
                epochStr = datestr(elements.epoch, 'yyyy-mm-ddTHH:MM:SS.FFF');
            else
                epochStr = elements.epoch;
            end
            if ~setXMLNodeValue(xmlDoc, 'EPOCH', epochStr, 'string')
                warning('EPOCH tag not found in XML file');
            end
        end

        % Update the six required orbital elements
        if ~setXMLNodeValue(xmlDoc, 'MEAN_MOTION', elements.meanMotion, 'numeric')
            warning('MEAN_MOTION tag not found in XML file');
        end
        if ~setXMLNodeValue(xmlDoc, 'ECCENTRICITY', elements.eccentricity, 'numeric')
            warning('ECCENTRICITY tag not found in XML file');
        end
        if ~setXMLNodeValue(xmlDoc, 'INCLINATION', elements.inclination, 'numeric')
            warning('INCLINATION tag not found in XML file');
        end
        if ~setXMLNodeValue(xmlDoc, 'RA_OF_ASC_NODE', elements.raan, 'numeric')
            warning('RA_OF_ASC_NODE tag not found in XML file');
        end
        if ~setXMLNodeValue(xmlDoc, 'ARG_OF_PERICENTER', elements.argPeriapsis, 'numeric')
            warning('ARG_OF_PERICENTER tag not found in XML file');
        end
        if ~setXMLNodeValue(xmlDoc, 'MEAN_ANOMALY', elements.meanAnomaly, 'numeric')
            warning('MEAN_ANOMALY tag not found in XML file');
        end

        userNodes = xmlDoc.getElementsByTagName('USER_DEFINED');
        for k = 0:userNodes.getLength-1
            node = userNodes.item(k);
            paramName = char(node.getAttribute('parameter'));
            if strcmp(paramName, 'SEMIMAJOR_AXIS')
                updateUserDefinedNode(node, a_km);
            elseif strcmp(paramName, 'PERIOD')
                updateUserDefinedNode(node, period_min);
            elseif strcmp(paramName, 'APOAPSIS')
                updateUserDefinedNode(node, apoapsis_alt);
            elseif strcmp(paramName, 'PERIAPSIS')
                updateUserDefinedNode(node, periapsis_alt);
            end
        end

        % Save the updated XML
        saveXML(xmlDoc, outputFile);

    catch ME
        error('Failed to update OMM file: %s', ME.message);
    end
end

function success = setXMLNodeValue(xmlDoc, tagName, newValue, valueType)
    success = false;
    nodes = xmlDoc.getElementsByTagName(tagName);
    if nodes.getLength() > 0
        node = nodes.item(0);
        
        % Format the value based on type
        if strcmp(valueType, 'string')
            valueStr = char(newValue);
        else
            valueStr = num2str(newValue, '%.8f');
        end
        
        if ~isempty(node.getFirstChild)
            node.getFirstChild.setNodeValue(valueStr);
            success = true;
        else
            % Create text node if it doesn't exist
            textNode = xmlDoc.createTextNode(valueStr);
            node.appendChild(textNode);
            success = true;
        end
    end
end

function updateUserDefinedNode(node, newValue)
    valueStr = num2str(newValue, '%.6f');
    if ~isempty(node.getFirstChild)
        node.getFirstChild.setNodeValue(valueStr);
    else
        textNode = node.getOwnerDocument.createTextNode(valueStr);
        node.appendChild(textNode);
    end
end

function saveXML(xmlDoc, filename)
    try
        import javax.xml.transform.*
        import javax.xml.transform.dom.*
        import javax.xml.transform.stream.*
        
        transformer = TransformerFactory.newInstance.newTransformer;
        transformer.setOutputProperty(OutputKeys.INDENT,'yes');
        transformer.setOutputProperty('{http://xml.apache.org/xslt}indent-amount','2');
        transformer.setOutputProperty(OutputKeys.OMIT_XML_DECLARATION,'no');
        
        source = DOMSource(xmlDoc);
        result = StreamResult(java.io.File(filename));
        transformer.transform(source,result);
    catch ME
        error('Failed to save XML file: %s', ME.message);
    end
end